CC = cc 
CCFLAGS = 
LDLIBS =  -lc -lm
OBJS = brusselator.o ranmar.o
SRCS = brusselator.c ranmar.c
EXE = simu
RM = /bin/rm -f

$(EXE): $(OBJS)
	$(CC) $(CCFLAGS) $(OBJS) -o $(EXE) $(LDLIBS)

.c.o:
	$(CC) $(CCFLAGS) -c $<

clean:
	$(RM) $(OBJS)

depend:
	makedepend $(SRCS)
