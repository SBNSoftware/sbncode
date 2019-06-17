# CAFAna starts all its includes with just CAFAna/ -- allow that
include_directories(..)

# And expects assert() to still function even in the prof build
add_definitions(-UNDEBUG)

add_subdirectory(StandardRecord)

# add_subdirectory(Core)
