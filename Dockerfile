FROM maven:3.9-eclipse-temurin-21 AS build
WORKDIR /app/ngspca
COPY ngspca/pom.xml .
COPY ngspca/src ./src
RUN mvn -B package \
 && cp target/ngspca-*.jar target/ngspca.jar

FROM eclipse-temurin:21-jre-jammy
WORKDIR /app

# --- Bundled tools for end-to-end 1000G high-coverage pipeline -----------

# mosdepth — fast BAM/CRAM depth calculation (static binary), pinned by
# checksum so a retagged release cannot change what the image computes
ARG MOSDEPTH_VERSION=0.3.14
ARG MOSDEPTH_SHA256=c5182b74a8f1b66710efa16e122cbc8a197834874b103e7c5c0bd9a6265ae7b6
ADD https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth \
    /usr/local/bin/mosdepth
RUN echo "${MOSDEPTH_SHA256}  /usr/local/bin/mosdepth" | sha256sum -c - \
 && chmod +x /usr/local/bin/mosdepth \
 && mosdepth --version

# --- NGS-PCA application ------------------------------------------------

COPY --from=build /app/ngspca/target/ngspca.jar /app/ngspca.jar
COPY resources /app/resources
ENV JAVA_TOOL_OPTIONS="-XX:+UseContainerSupport -XX:MaxRAMPercentage=90.0"
ENTRYPOINT ["java", "-jar", "/app/ngspca.jar"]
