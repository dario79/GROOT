FROM alexo91/groot-soil:latest as builder0
LABEL stage="builder"

USER root

#compile and install GROOT
RUN mkdir -p /opt/GROOT/source \
    && mkdir /opt/GROOT/build \
    && mkdir /opt/GROOT/install 

RUN git clone https://github.com/dario79/GROOT.git /opt/GROOT/source/
COPY ./build-GROOT.sh /opt/GROOT/source/

RUN chmod +x /opt/GROOT/source/build-GROOT.sh
RUN "/opt/GROOT/source/build-GROOT.sh"

#--------------------------------------------------

FROM alexo91/groot-soil:latest

COPY --from=builder0 --chown=groot:groot /opt/GROOT/install /opt/GROOT/
COPY --chown=groot:groot ./entrypoint.sh /


RUN chmod +x ./entrypoint.sh
ENTRYPOINT ["./entrypoint.sh"]
CMD ["GROOT"]