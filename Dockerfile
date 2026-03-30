FROM maven:3.9-eclipse-temurin-11 AS build
WORKDIR /app/ngspca
COPY ngspca/pom.xml .
COPY ngspca/src ./src
RUN mvn -B package \
 && cp target/ngspca-*.jar target/ngspca.jar

FROM eclipse-temurin:11-jre-jammy
WORKDIR /app

# --- Bundled tools for end-to-end 1000G high-coverage pipeline -----------

# mosdepth v0.3.9 — fast BAM/CRAM depth calculation (static binary)
ARG MOSDEPTH_VERSION=0.3.9
ADD https://github.com/brentp/mosdepth/releases/download/v${MOSDEPTH_VERSION}/mosdepth \
    /usr/local/bin/mosdepth
RUN chmod +x /usr/local/bin/mosdepth

# --- NGS-PCA application ------------------------------------------------

COPY --from=build /app/ngspca/target/ngspca.jar /app/ngspca.jar
COPY resources /app/resources
ENV JAVA_TOOL_OPTIONS="-XX:+UseContainerSupport -XX:MaxRAMPercentage=90.0"
ENTRYPOINT ["java", "-jar", "/app/ngspca.jar"]
