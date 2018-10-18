###################################################################
# Package: RJTupler 
###################################################################

# declare the package name
atlas_subdir( RJTupler )

# extra deps based on the build environment:
set( extra_deps )
set( extra_libs )

# declare the RJTupler package dependencies
atlas_depends_on_subdirs(
    PUBLIC
    Superflow
    SusyNtuple
    ${extra_deps}
)

#external dependencies
find_package( ROOT COMPONENTS Gpad Graf Graf3d Core Tree MathCore Hist RIO )

include_directories($ENV{TestArea}/../RestFrames/include)
find_library(RJLIB RestFrames)

#atlas_add_root_dictionary( RJTuplerLib RJTuplerCintDict
#        ROOT_HEADERS  Root/LinkDef.h
#        EXTERNAL_PACKAGES ROOT
#)

# libraries
atlas_add_library( RJTuplerLib
    RJTupler/*.h Root/*.cxx ${RJTuplerCintDict}
    PUBLIC_HEADERS RJTupler
    INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} #$ENV{TestArea}/../RestFrames/include
    LINK_LIBRARIES SuperflowLib SusyNtupleLib ${RJLIB} ${ROOT_LIBRARIES}
)

# executable(s) in the package
set( extra_libs )

# define a function to call add_executable macro on our executables
function( SusyNtExec filename)
    set(execname)
    get_filename_component(execname ${filename} NAME_WE)
    atlas_add_executable( ${execname} "util/${execname}.cxx"
        INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
        LINK_LIBRARIES ${ROOT_LIBRARIES} RJTuplerLib ${extra_libs}
    )
endfunction( SusyNtExec )

# build everything in util/ as an executable for RJTupler
file(GLOB files "util/*.cxx")
foreach(file ${files})
    SusyNtExec(${file})
endforeach()

atlas_install_data( data/* )
