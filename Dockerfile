FROM debian:12.4
RUN apt update && \
    apt install -y git cmake clang libboost-all-dev nlohmann-json3-dev
WORKDIR /src
COPY . .
WORKDIR /app
RUN cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ /src && \
    cmake --build /app --target ai_vp_cpp --config Release -v
EXPOSE 8080
CMD ["ai_vp_cpp/ai_vp_cpp"]
