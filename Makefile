PROG = rstat

SRCS = rstat.c
HEADERS =
OBJS = $(SRCS:.c=.o)

CC = gcc
yacc = yacc
CFLAGS = -Wall -g -O2
LDFLAGS += -lm

all: $(PROG)

clean:
	rm -f $(OBJS) $(PROG)

$(PROG): $(OBJS) $(HEADERS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)
