FROM ubuntu:16.04
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates curl gcc gcc-arm-none-eabi libc6-dev libcurl4-openssl-dev libssh2-1 libnewlib-dev qemu-user-static
RUN curl -LSfs https://japaric.github.io/trust/install.sh | \
    sh -s -- --git japaric/xargo --tag v0.3.1 --target x86_64-unknown-linux-gnu --to /usr/bin
ENV AR_thumbv6m_linux_eabi=arm-none-eabi-ar \
    CARGO_TARGET_THUMBV6M_LINUX_EABI_LINKER=arm-none-eabi-gcc \
    CC_thumbv6m_linux_eabi=arm-none-eabi-gcc \
