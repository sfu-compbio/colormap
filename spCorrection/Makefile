CC=g++

OEACORR_VERSION := "0.1"
BUILD_DATE := "$(shell date)"

all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS build
build: clean spCorrection

LDFLAGS=#-static
LIBS=-lz -lm -pthread -lpthread
CFLAGS=-DSPCORR_VERSION=\"$(SPCORR_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" #-fno-pic 

objects=main.o common.o reader.o correct.o

# compile: $(objects)

spCorrection: $(objects)
	$(CC) $(objects) -o $@ ${LDFLAGS} ${LIBS}

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	@rm -f $(objects) spCorrection
	
clean-executable:
	@rm -f spCorrection

DEBUG_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -ggdb)
	$(eval LIBS = $(LIBS) -ggdb)

OPTIMIZE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -O2)

PROFILE_FLAGS:
		$(eval CFLAGS = $(CFLAGS) -pg -g)
	$(eval LIBS = $(LIBS) -pg -g)
