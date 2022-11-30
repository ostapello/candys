QT += widgets opengl
CONFIG += c++1z

SOURCES += \
    buttonsarea.cpp \
    functions.cpp \
    grapharea.cpp \
    main.cpp \
    window.cpp

HEADERS += \
    buttonsarea.h \
    functions.h \
    grapharea.h \
    window.h

QMAKE_CXXFLAGS += -Wextra -pedantic -Werror
