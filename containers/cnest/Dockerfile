FROM tomas81/cnest@sha256:092d680fcc737c92056abf16aa3e344a859b83aa296020e30c0da6afbab6cb5c

COPY CNest/src/cnest_dev.py /resources/cnest_dev.py
RUN chmod +x /resources/cnest_dev.py

# Unset entrypoint
ENTRYPOINT []

# Set default command
CMD ["bash"]
