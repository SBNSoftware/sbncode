#
# NOTE (from Gray Putnam):
#
# In order to get sbn to work on the grid, we need dictionaries to be defined with a relative 
# path to the source header file. The current ROOT_GENERATE_DICTIONARY CMake macro defines
# an absolute path.
#
# This behavior is changed in ROOT v6.18 to use a relative path (which we want). So this
# is a copy of ROOT_GENERATE_DICTIONARY on that version. Use it until we upgrade by replacing
# "ROOT_GENERATE_DICTIONARY" with "FUTURE_ROOT_GENERATE_DICTIONARY"


#---------------------------------------------------------------------------------------------------
#---FUTURE_ROOT_GENERATE_DICTIONARY( result_var )
# Returns the path to the .so file or .dll file. In the latter case Windows defines the dll files as
# executables and puts them in the $ROOTSYS/bin folder.
function(ROOT_GET_LIBRARY_OUTPUT_DIR result)
  set(library_output_dir)
  if(MSVC)
    if(DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY AND NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY STREQUAL "")
      set(library_output_dir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    else()
      set(library_output_dir ${CMAKE_CURRENT_BINARY_DIR})
    endif()
  else()
    if(DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY AND NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY STREQUAL "")
      set(library_output_dir ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    else()
      set(library_output_dir ${CMAKE_CURRENT_BINARY_DIR})
    endif()
  endif()
  SET(${result} "${library_output_dir}" PARENT_SCOPE)
endfunction(ROOT_GET_LIBRARY_OUTPUT_DIR)
#---------------------------------------------------------------------------------------------------
#---FUTURE_ROOT_GENERATE_DICTIONARY( dictionary headerfiles NODEPHEADERS ghdr1 ghdr2 ...
#                                                    MODULE module DEPENDENCIES dep1 dep2
#                                                    BUILTINS dep1 dep2
#                                                    STAGE1 LINKDEF linkdef OPTIONS opt1 opt2 ...)
#
# <dictionary> is the dictionary stem; the macro creates (among other files) the dictionary source as
#   <dictionary>.cxx
# <headerfiles> are "as included"; set appropriate INCLUDE_DIRECTORIES property on the directory.
#   The dictionary target depends on these headers. These files must exist.
# <NODEPHEADERS> same as <headerfiles>. If these files are not found (given the target include path)
#   no error is emitted. The dictionary does not depend on these headers.
#---------------------------------------------------------------------------------------------------
function(FUTURE_ROOT_GENERATE_DICTIONARY dictionary)
  CMAKE_PARSE_ARGUMENTS(ARG "STAGE1;MULTIDICT;NOINSTALL" "MODULE;LINKDEF" "NODEPHEADERS;OPTIONS;DEPENDENCIES;BUILTINS" ${ARGN})

  # Check if OPTIONS start with a dash.
  if (ARG_OPTIONS)
    foreach(ARG_O ${ARG_OPTIONS})
      if (NOT ARG_O MATCHES "^-*")
        message(FATAL_ERROR "Wrong rootcling option: ${ARG_OPTIONS}")
      endif()
    endforeach()
  endif(ARG_OPTIONS)

  #---roottest compability---------------------------------
  if(CMAKE_ROOTTEST_DICT)
    set(CMAKE_INSTALL_LIBDIR ${CMAKE_CURRENT_BINARY_DIR})
    set(libprefix "")
  endif()

  #---Get the list of include directories------------------
  get_directory_property(incdirs INCLUDE_DIRECTORIES)
  # rootcling invoked on foo.h should find foo.h in the current source dir,
  # no matter what.
  list(APPEND incdirs ${CMAKE_CURRENT_SOURCE_DIR})

  #---Get the list of header files-------------------------
  # CMake needs dependencies from ${CMAKE_CURRENT_SOURCE_DIR} while rootcling wants
  # header files "as included" (and thus as passed as argument to this CMake function).
  set(headerfiles)
  set(_list_of_header_dependencies)
  foreach(fp ${ARG_UNPARSED_ARGUMENTS})
    if(${fp} MATCHES "[*?]") # Is this header a globbing expression?
      file(GLOB files inc/${fp} ${fp}) # Elements of ${fp} have the complete path.
      foreach(f ${files})
        if(NOT f MATCHES LinkDef) # skip LinkDefs from globbing result
          set(add_inc_as_include On)
          string(REGEX REPLACE "^${CMAKE_CURRENT_SOURCE_DIR}/inc/" "" f_no_inc ${f})
          list(APPEND headerfiles ${f_no_inc})
          list(APPEND _list_of_header_dependencies ${f})
        endif()
      endforeach()
    else()
      if(IS_ABSOLUTE ${fp})
        set(headerFile ${fp})
      else()
        find_file(headerFile ${fp} HINTS ${incdirs} NO_DEFAULT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
      endif()
      if(NOT headerFile)
        message(FATAL_ERROR "Cannot find header ${fp} to generate dictionary ${dictionary} for. Did you forget to set the INCLUDE_DIRECTORIES property for the current directory?")
      endif()
      list(APPEND headerfiles ${fp})
      list(APPEND _list_of_header_dependencies ${headerFile})
      unset(headerFile CACHE) # find_file, forget headerFile!
    endif()
  endforeach()

  foreach(fp ${ARG_NODEPHEADERS})
    list(APPEND headerfiles ${fp})
    # no dependency - think "vector" etc.
  endforeach()

  if(NOT (headerfiles OR ARG_LINKDEF))
    message(FATAL_ERROR "No headers nor LinkDef.h supplied / found for dictionary ${dictionary}!")
  endif()

  if(CMAKE_PROJECT_NAME STREQUAL ROOT)
    set(includedirs -I${CMAKE_SOURCE_DIR}
                    -I${CMAKE_BINARY_DIR}/etc/cling/ # This is for the RuntimeUniverse
                    -I${CMAKE_BINARY_DIR}/include)
    set(excludepaths ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/inc)
    set(includedirs -I${CMAKE_CURRENT_SOURCE_DIR}/inc)
  endif()
  foreach( d ${incdirs})
   set(includedirs ${includedirs} -I${d})
  endforeach()

  foreach(dep ${ARG_DEPENDENCIES})
    if(TARGET ${dep})
      get_property(dep_include_dirs TARGET ${dep} PROPERTY INCLUDE_DIRECTORIES)
      foreach(d ${dep_include_dirs})
        set(includedirs ${includedirs} -I${d})
      endforeach()
    endif()
  endforeach()

  if(includedirs)
    list(REMOVE_DUPLICATES includedirs)
  endif()

  #---Get the list of definitions---------------------------
  get_directory_property(defs COMPILE_DEFINITIONS)
  foreach( d ${defs})
   if((NOT d MATCHES "=") AND (NOT d MATCHES "^[$]<.*>$")) # avoid generator expressions
     set(definitions ${definitions} -D${d})
   endif()
  endforeach()
  #---Get LinkDef.h file------------------------------------
  foreach( f ${ARG_LINKDEF})
    if( IS_ABSOLUTE ${f})
      set(_linkdef ${_linkdef} ${f})
    else()
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/inc/${f})
        set(_linkdef ${_linkdef} ${CMAKE_CURRENT_SOURCE_DIR}/inc/${f})
      else()
        set(_linkdef ${_linkdef} ${CMAKE_CURRENT_SOURCE_DIR}/${f})
      endif()
    endif()
  endforeach()

  #---Build the names for library, pcm and rootmap file ----
  get_filename_component(dict_base_name ${dictionary} NAME_WE)
  if(dict_base_name MATCHES "^G__")
    string(SUBSTRING ${dictionary} 3 -1 deduced_arg_module)
  else()
    set(deduced_arg_module ${dict_base_name})
  endif()

  #---Set the library output directory-----------------------
  ROOT_GET_LIBRARY_OUTPUT_DIR(library_output_dir)
  if(ARG_MODULE)
    set(library_name ${libprefix}${ARG_MODULE}${libsuffix})
    if(ARG_MULTIDICT)
      set(newargs -s ${library_output_dir}/${library_name} -multiDict)
      set(pcm_name ${library_output_dir}/${libprefix}${ARG_MODULE}_${dictionary}_rdict.pcm)
      set(rootmap_name ${library_output_dir}/${libprefix}${deduced_arg_module}.rootmap)
    else()
      set(newargs -s ${library_output_dir}/${library_name})
      set(pcm_name ${library_output_dir}/${libprefix}${ARG_MODULE}_rdict.pcm)
      set(rootmap_name ${library_output_dir}/${libprefix}${ARG_MODULE}.rootmap)
      set(cpp_module ${ARG_MODULE})
      if(runtime_cxxmodules)
        set(cpp_module_file ${library_output_dir}/${cpp_module}.pcm)
        if (APPLE)
          if (${cpp_module} MATCHES "(GCocoa|GQuartz)")
            set(cpp_module_file)
          endif()
        endif(APPLE)
      endif()
    endif()
  else()
    set(library_name ${libprefix}${deduced_arg_module}${libsuffix})
    set(newargs -s ${library_output_dir}/${library_name})
    set(pcm_name ${library_output_dir}/${libprefix}${deduced_arg_module}_rdict.pcm)
    set(rootmap_name ${library_output_dir}/${libprefix}${deduced_arg_module}.rootmap)
  endif()

  if(CMAKE_ROOTTEST_NOROOTMAP)
    set(rootmap_name )
    set(rootmapargs )
  else()
    set(rootmapargs -rml ${library_name} -rmf ${rootmap_name})
  endif()

  #---Get the library and module dependencies-----------------
  if(ARG_DEPENDENCIES)
    foreach(dep ${ARG_DEPENDENCIES})
      set(newargs ${newargs} -m  ${libprefix}${dep}_rdict.pcm)
    endforeach()
  endif()

  if(cpp_module_file)
    set(newargs -cxxmodule ${newargs})
  endif()

  #---what rootcling command to use--------------------------
  if(ARG_STAGE1)
    set(command ${CMAKE_COMMAND} -E env "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}" $<TARGET_FILE:rootcling_stage1>)
    set(pcm_name)
  else()
    if(CMAKE_PROJECT_NAME STREQUAL ROOT)
      set(command ${CMAKE_COMMAND} -E env "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib:$ENV{LD_LIBRARY_PATH}"
                  "ROOTIGNOREPREFIX=1" $<TARGET_FILE:rootcling> -rootbuild)
      set(ROOTCINTDEP rootcling)
    elseif(TARGET ROOT::rootcling)
      set(command ${CMAKE_COMMAND} -E env "LD_LIBRARY_PATH=${ROOT_LIBRARY_DIR}:$ENV{LD_LIBRARY_PATH}" $<TARGET_FILE:ROOT::rootcling>)
    else()
      set(command ${CMAKE_COMMAND} -E env rootcling)
    endif()
  endif()

  #---build the path exclusion switches----------------------
  set(excludepathsargs "")
  foreach(excludepath ${excludepaths})
    set(excludepathsargs ${excludepathsargs} -excludePath ${excludepath})
  endforeach()

  #---build the implicit dependencies arguments
  # NOTE: only the Makefile generator respects this!
  foreach(_dep ${_linkdef} ${_list_of_header_dependencies})
    list(APPEND _implicitdeps CXX ${_dep})
  endforeach()

  if(ARG_MODULE)
    set(MODULE_LIB_DEPENDENCY ${ARG_DEPENDENCIES})
  endif()

  #---call rootcint------------------------------------------
  add_custom_command(OUTPUT ${dictionary}.cxx ${pcm_name} ${rootmap_name} ${cpp_module_file}
                     COMMAND ${command} -v2 -f  ${dictionary}.cxx ${newargs} ${excludepathsargs} ${rootmapargs}
                                        ${definitions} ${includedirs} ${ARG_OPTIONS} ${headerfiles} ${_linkdef}
                     IMPLICIT_DEPENDS ${_implicitdeps}
                     DEPENDS ${_list_of_header_dependencies} ${_linkdef} ${ROOTCINTDEP} ${MODULE_LIB_DEPENDENCY})
  get_filename_component(dictname ${dictionary} NAME)

  #---roottest compability
  add_custom_target(${dictname} DEPENDS ${dictionary}.cxx ${pcm_name} ${rootmap_name} ${cpp_module_file})

  if(NOT ARG_NOINSTALL AND NOT CMAKE_ROOTTEST_DICT AND DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set_property(GLOBAL APPEND PROPERTY ROOT_DICTIONARY_TARGETS ${dictname})
    set_property(GLOBAL APPEND PROPERTY ROOT_DICTIONARY_FILES ${CMAKE_CURRENT_BINARY_DIR}/${dictionary}.cxx)

    ROOT_GET_INSTALL_DIR(shared_lib_install_dir)
    # Install the C++ module if we generated one.
    if (cpp_module_file)
      install(FILES ${cpp_module_file}
                    DESTINATION ${shared_lib_install_dir} COMPONENT libraries)
    endif()

    if(ARG_STAGE1)
      install(FILES ${rootmap_name}
                    DESTINATION ${shared_lib_install_dir} COMPONENT libraries)
    else()
      install(FILES ${pcm_name} ${rootmap_name}
                    DESTINATION ${shared_lib_install_dir} COMPONENT libraries)
    endif()
  endif()

  if(ARG_BUILTINS)
    foreach(arg1 ${ARG_BUILTINS})
      if(${arg1}_TARGET)
        add_dependencies(${dictname} ${${arg1}_TARGET})
      endif()
    endforeach()
  endif()
  # FIXME: Support mulptiple dictionaries. In some cases (libSMatrix and
  # libGenVector) we have to have two or more dictionaries (eg. for math,
  # we need the two for double vs Double32_t template specializations).
  # In some other cases, eg. libTreePlayer.so we add in a separate dictionary
  # files which for some reason (temporarily?) cannot be put in the PCH. Eg.
  # all rest of the first dict is in the PCH but this file is not and it
  # cannot be present in the original dictionary.
  if(cpp_module)
    ROOT_CXXMODULES_APPEND_TO_MODULEMAP("${cpp_module}" "${headerfiles}")
  endif()
endfunction(FUTURE_ROOT_GENERATE_DICTIONARY)

