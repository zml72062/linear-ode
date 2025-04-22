SRCDIR    = src
INCDIR    = include
OBJDIR    = build
CXX       = g++
CPPFLAGS  = -I./${INCDIR}
CXXFLAGS  = -O2 -Wall
LDFLAGS   = -lcln -lginac -lflint


OBJS      = ${OBJDIR}/iszero.cpp.o \
			${OBJDIR}/interface.cpp.o \
			${OBJDIR}/base.cpp.o \
			${OBJDIR}/diffeq.cpp.o \
			${OBJDIR}/solver.cpp.o

all: pre example

example: ${OBJS} example.cpp
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${OBJS} example.cpp -o example ${LDFLAGS}

${OBJS}: ${OBJDIR}/%.cpp.o: ${SRCDIR}/%.cpp
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $^ -o $@

.PHONY: pre
pre:
	mkdir -p ${OBJDIR}

.PHONY: clean
clean:
	rm -rf ${OBJDIR}
	rm -f example


