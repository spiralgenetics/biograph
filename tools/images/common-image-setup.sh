#!/bin/bash
set -e

# Common commands to be executed to set up a new image.  This is run
# before any image-specific scripts are run.

export RELEASE=`/usr/bin/lsb_release -sc`

# Avoid prompts from everything
export DEBIAN_FRONTEND=noninteractive

cat <<EOF > /etc/rsyslog.d/S99-spiral.conf
# Spiral services

local0.info		/var/log/spiral
local0.debug	/var/log/spiral.debug
EOF

echo "search spiralgenetics.com" > /etc/resolvconf/resolv.conf.d/base
echo "kernel.core_pattern = /tmp/core.%h.%e.%p" > /etc/sysctl.d/50-core.conf
service procps restart

# echo 'APT::Get::Install-Recommends "false";' > /etc/apt/apt.conf

apt-get -yq update
apt-get -yq upgrade
apt-get -yq install \
    bash-completion \
    bc \
    curl \
    htop \
    iotop \
    iptraf-ng \
    fio \
    jq \
    less \
    libbz2-dev \
    libcurl3 \
    libdaemon0 \
    libffi-dev \
    liblzma-dev \
    libgomp1 \
    libssl-dev \
    libxml2 \
    lnav \
    ipython \
    ncbi-blast+ \
    ncbi-data \
    nginx \
    nfs-kernel-server \
    python-meld3 \
    python-virtualenv \
    python-dev \
    python-pip \
    libffi-dev \
    libssl-dev \
    samtools \
    supervisor \
    tree \
    unzip \
    vim

# Not sure where this comes from in Trusty 14.04.5 but it's in the way.
rm -f /usr/lib/python*/dist-packages/pyOpenSSL-0.13.egg-info

# Specific version of awscli needed because syntax changes between releases
for pkg in setuptools Markdown pyOpenSSL cryptography pyvcf networkx pybedtools 'awscli==1.11.28' swalign 'pylint==1.7.2' progressbar2; do
  pip install --upgrade $pkg
done

cat <<EOF >/etc/logrotate.d/spiral
/var/log/spiral {
	rotate 14
	daily
	compress
	missingok
	notifempty
	postrotate
		reload rsyslog >/dev/null 2>&1 || true
	endscript
}

/var/log/spiral.debug {
	rotate 14
	daily
	compress
	missingok
	notifempty
	postrotate
		reload rsyslog >/dev/null 2>&1 || true
	endscript
}
EOF
