FROM ubuntu:18.04
 
# Updates
RUN apt-get update && apt-get install -y cmake libboost-all-dev g++

# Make up to build
RUN mkdir -p /root/htstream/build

# Copy current directory (./HTStream/.) to the correct place.
# If we used CI/CD - it will be programmatically set to master or the proper branch to test.
# If it is local, it _might_ be off, but it is nice for developers. However, it is probably
# nicer for developers if we just mount these in. But it might be more robust if we just download git
# and clone it. Happy to change this if someone has really strong feelings, I do not.
COPY . /root/htstream/

# We do our build in `./htstream/build` so it makes sense to put it here.
WORKDIR /root/htstream/build

# Build HTStream all in the same layer. 
RUN cmake .. && \
	make && \
	make install
