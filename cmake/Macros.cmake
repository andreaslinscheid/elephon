
macro( prefer_static )
    if( NOT WIN32 )
        list( REMOVE_ITEM CMAKE_FIND_LIBRARY_SUFFIXES   ".a" )
        list( INSERT      CMAKE_FIND_LIBRARY_SUFFIXES 0 ".a" )
    endif()
endmacro()

macro( prefer_dynamic )
    if( NOT WIN32 )
        list( REMOVE_ITEM CMAKE_FIND_LIBRARY_SUFFIXES ".a" )
        list( APPEND      CMAKE_FIND_LIBRARY_SUFFIXES ".a" )
    endif()
endmacro()
