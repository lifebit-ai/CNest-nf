FROM quay.io/lifebitaiorg/cnest_tomas81:dev2023_april

COPY CNest/src/cnest_dev.py /resources/cnest_dev.py
RUN chmod +x /resources/cnest_dev.py

# Unset entrypoint
ENTRYPOINT []

# Set default command
CMD ["bash"]
