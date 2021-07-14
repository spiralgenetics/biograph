FROM python:3.7
WORKDIR /

ADD app /app
ADD lib /mylib

# Install any necessary dependencies
RUN apt-get update
RUN apt-get install -y bcftools tabix samtools

RUN pip install /app/BioGraph*.tgz

ENV LD_LIBRARY_PATH /mylib
ENV PATH /app:/bin:/usr/bin:/usr/local/bin

CMD ["cat", "/app/biograph-help.txt"]
