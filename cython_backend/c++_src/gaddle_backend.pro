#-------------------------------------------------
#
# Project created by QtCreator 2017-11-04T20:57:11
#
#-------------------------------------------------

TARGET = gaddle_backend
TEMPLATE = lib

DEFINES += GADDLE_BACKEND_LIBRARY

SOURCES += gaddle_backend.cpp \
    algorithm_utils.cpp \
    transform_molecule.cpp \
    minimize.cpp \

HEADERS += gaddle_backend.h\
        gaddle_backend_global.h \
    algorithm_utils.h \
    common_libraries.h \
    transform_molecule.h \
    minimize.h

CONFIG += c++11

unix {
    target.path = /usr/lib
    INSTALLS += target
}

unix|win32: LIBS += -larmadillo
