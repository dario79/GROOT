FROM alexo91/alma9g4r-base:latest as builder0
LABEL stage="builder"


RUN mkdir -p /opt/GROOT/source \
    && mkdir /opt/GROOT/build \
    && mkdir /opt/GROOT/install 

# COPY ./GROOT.cc /opt/GROOT/source/
# COPY ./CMakeLists.txt /opt/GROOT/source/
# COPY ./inputFileDetectors.txt /opt/GROOT/source/
# COPY ./run.mac /opt/GROOT/source/
# COPY ./vis.mac /opt/GROOT/source/
# COPY ./src/ /opt/GROOT/source/src/
# COPY ./include/ /opt/GROOT/source/include/
# COPY ./cmake/ /opt/GROOT/source/cmake/

RUN git clone https://github.com/dario79/GROOT.git /opt/GROOT/source/
COPY ./build-GROOT.sh /opt/GROOT/source/

RUN chmod +x /opt/GROOT/source/build-GROOT.sh
RUN "/opt/GROOT/source/build-GROOT.sh"

#--------------------------------------------------

FROM alexo91/alma9g4r-base:latest

COPY --from=builder0 /opt/GROOT/install /opt/GROOT/

COPY ./entrypoint.sh /
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
