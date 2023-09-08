PROG = rstat

SRCS = rstat.cc
HEADERS =
OBJS = $(SRCS:.c=.o)

CC=clang++
CFLAGS = -Wall -g
LDFLAGS += -lm

all: $(PROG)

clean:
	rm -f $(OBJS) $(PROG)

$(PROG): $(OBJS) $(HEADERS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)
