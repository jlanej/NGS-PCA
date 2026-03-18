FROM maven:3.6.3-openjdk-11 AS build
WORKDIR /app/ngspca
COPY ngspca/pom.xml .
COPY ngspca/src ./src
RUN mvn install

FROM maven:3.6.3-openjdk-11
WORKDIR /app

# --- Bundled tools for end-to-end 1000G high-coverage pipeline -----------

# mosdepth v0.3.9 — fast BAM/CRAM depth calculation (static binary)
ARG MOSDEPTH_VERSION=0.3.9
ADD https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth \
    /usr/local/bin/mosdepth
RUN chmod +x /usr/local/bin/mosdepth

# IBM Aspera Connect — high-speed download from EBI/NCBI
ARG ASPERA_VERSION=4.2.12.780
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends wget && \
    wget -qO /tmp/aspera.tar.gz \
      "https://d3gcli72yxqn2z.cloudfront.net/downloads/connect/latest/bin/ibm-aspera-connect_${ASPERA_VERSION}_linux_x86_64.tar.gz" && \
    tar xzf /tmp/aspera.tar.gz -C /tmp && \
    bash /tmp/ibm-aspera-connect_${ASPERA_VERSION}_linux_x86_64.sh && \
    ln -sf /root/.aspera/connect/bin/ascp /usr/local/bin/ascp && \
    rm -rf /tmp/* /var/lib/apt/lists/*

ENV ASPERA_SSH_KEY=/root/.aspera/connect/etc/asperaweb_id_dsa.openssh

# --- NGS-PCA application ------------------------------------------------

COPY --from=build /app/ngspca/target/ngspca-0.02-SNAPSHOT.jar /app
COPY resources /app/resources
ENV JAVA_TOOL_OPTIONS "-XX:+UseContainerSupport -XX:MaxRAMPercentage=90.0"
ENTRYPOINT ["java", "-jar","ngspca-0.02-SNAPSHOT.jar"]
