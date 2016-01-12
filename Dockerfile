FROM kbase/kbase:sdkbase.latest
MAINTAINER Fangfang Xia
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.

# RUN apt-get update

# -----------------------------------------


RUN apt-get install libffi-dev libssl-dev
RUN pip install --upgrade requests[security]

# Install AssemblyRAST client

RUN pip install PrettyTable

RUN \
    git clone https://github.com/kbase/assembly.git && \
    cd assembly && \
    git checkout next && \
    make -f Makefile.standalone


# Copy local wrapper files, and build

COPY ./ /kb/module
RUN mkdir -p /kb/module/work

WORKDIR /kb/module

RUN make

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
