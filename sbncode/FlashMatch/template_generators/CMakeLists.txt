# All of the .py files in this directory are python modules.
# We distinguish executable and non-executable .py files.
# The former are installed with execute permission.
# The latter are installed without execute permission.

# Executable python files.

LIST(APPEND exes generate_simple_weighted_template.py )

# Non-executable python files.

# LIST(APPEND nonexes experiment_utilities.py )

message(STATUS "Executable python modules ${exes}")
# message(STATUS "Non-executable python modules ${nonexes}")

install(PROGRAMS ${exes} DESTINATION ${product}/${version}/python)
# install(FILES ${nonexes} DESTINATION ${product}/${version}/python)
install_scripts(LIST ${exes})
