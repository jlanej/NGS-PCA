FROM maven:3.6.3-openjdk-11 AS build
WORKDIR /app/ngspca
COPY ngspca/pom.xml .
COPY ngspca/src ./src
RUN mvn install

FROM maven:3.6.3-openjdk-11
WORKDIR /app
COPY --from=build /app/ngspca/target/ngspca-0.02-SNAPSHOT.jar /app
ENV JAVA_TOOL_OPTIONS "-XX:+UseContainerSupport -XX:MaxRAMPercentage=90.0"
ENTRYPOINT ["java", "-jar","ngspca-0.02-SNAPSHOT.jar"]
