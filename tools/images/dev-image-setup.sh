#!/bin/bash
export RELEASE=`/usr/bin/lsb_release -sc`

# Azure or AWS?
# http://blog.mszcool.com/index.php/2015/04/detecting-if-a-virtual-machine-runs-in-microsoft-azure-linux-windows-to-protect-your-software-when-distributed-via-the-azure-marketplace/
if `grep -q unknown-245 /var/lib/dhcp/dhclient.eth0.leases`; then
    export PLATFORM="Azure"
else
    export PLATFORM="AWS"
fi

# Install boost packages we depend on
if test -d /usr/local/include/boost
then
    echo "Unexpected boost installation present in /usr/local/include/"
    exit 1
fi

# Add bazel package repository (jdk 1.8 required)
echo "deb [arch=amd64] http://storage.googleapis.com/bazel-apt stable jdk1.8" > /etc/apt/sources.list.d/bazel.list

# Add bazel key
apt-key add - <<EOF
-----BEGIN PGP PUBLIC KEY BLOCK-----
Version: GnuPG v1

mQINBFdEmzkBEACzj8tMYUau9oFZWNDytcQWazEO6LrTTtdQ98d3JcnVyrpT16yg
I/QfGXA8LuDdKYpUDNjehLtBL3IZp4xe375Jh8v2IA2iQ5RXGN+lgKJ6rNwm15Kr
qYeCZlU9uQVpZuhKLXsWK6PleyQHjslNUN/HtykIlmMz4Nnl3orT7lMI5rsGCmk0
1Kth0DFh8SD9Vn2G4huddwxM8/tYj1QmWPCTgybATNuZ0L60INH8v6+J2jJzViVc
NRnR7mpouGmRy/rcr6eY9QieOwDou116TrVRFfcBRhocCI5b6uCRuhaqZ6Qs28Bx
4t5JVksXJ7fJoTy2B2s/rPx/8j4MDVEdU8b686ZDHbKYjaYBYEfBqePXScp8ndul
XWwS2lcedPihOUl6oQQYy59inWIpxi0agm0MXJAF1Bc3ToSQdHw/p0Y21kYxE2pg
EaUeElVccec5poAaHSPprUeej9bD9oIC4sMCsLs7eCQx2iP+cR7CItz6GQtuZrvS
PnKju1SKl5iwzfDQGpi6u6UAMFmc53EaH05naYDAigCueZ+/2rIaY358bECK6/VR
kyrBqpeq6VkWUeOkt03VqoPzrw4gEzRvfRtLj+D2j/pZCH3vyMYHzbaaXBv6AT0e
RmgtGo9I9BYqKSWlGEF0D+CQ3uZfOyovvrbYqNaHynFBtrx/ZkM82gMA5QARAQAB
tEdCYXplbCBEZXZlbG9wZXIgKEJhemVsIEFQVCByZXBvc2l0b3J5IGtleSkgPGJh
emVsLWRldkBnb29nbGVncm91cHMuY29tPokCPgQTAQIAKAUCV0SbOQIbAwUJA8Jn
AAYLCQgHAwIGFQgCCQoLBBYCAwECHgECF4AACgkQPVkZtEhFfuA+BA/9ErVPkCbr
bwGX6HHY2h9RULnOuBOJKhqDdlAmJ9Bi1TwliG72rXO52mO7WbSJ4Ux+5bbO3kiz
UId6xr6FZmTJiHCY42gwgLOHLIwd6Y7l+2PTyQ0ryL35hiPvlwhcJ4QzfF0/8MTD
vyQaieBKQQkZD+RmXJLmDGmztLhaUw0BTkONfbN0jTdIt7lAmkVZYo8JCoMbCXlk
9kkh6lx9vXZ0JQkfm1gRRWAAShsOHMYZgXsj11Y3C76zxv8jNoT8315Q4bVMKmjr
+bOaR7idiTlxRKl17eRr8qZ1AqJrwywyjFTWQpwHw3uR4iMt+QSgceKiuihpZBpv
mVtNVaMiJTkR63pRmwaGoGov4u3VMBRa6YMQaXQ1UEwrqEH+cyNZgWidlIywpgzD
4hjLYuE/dCvCVT+M32vTQVhNDCwZkmiNeVQE1R/5AdhQLSvBnLGXIfkXnsCl0vFl
oG2brgxIs3TN8oEitIeVWJ+oj5AqTRBeHIZ5MHktF/18AfPNER4UtufGxUWbAtyf
3HtmZZ03XMct52qSaKZbTY1bOmsgy6ba3Kkc9ifjYa9w9FEVysUPmLzy4+l3K5Xw
92iVYLDtCTThhYfbTJ8IWvJGE9hkyYmPF4Q0O3kFQqQmrFVADMyxR0DMIULs8U3k
bKa5K3939aSYXQAdpq6d6T0y2+XdvSHrmVu5Ag0EV0SbOQEQAOef9VQZQ6VfxJVM
i5kcjws/1fprB3Yp8sODL+QyULqbmcJTMr8Tz83OxprCH5Nc7jsw1oqzbNtq+N2p
OnbAL6XFPolQYuOjKlHGzbQvpH8ZSok6AzwrPNq3XwoB0+12A86wlpajUPfvgajN
jmESMchLnIs3qH1j5ayVICr7vH1i1Wem2J+C/z6gIaG4bko0XKAeU6fNYRmuHLHC
iBiKocpn54LmmPL4ifN7Rz1KkCaAKTT8vKtaVh0g1eswb+9W3qldm+nAc6e1ajWD
iLqhOmTQRVrght80XPYmtv2x8cdkxgECbT6T84rZtMZAdxhjdOmJ50ghPn9o/uxd
CDurhZUsu4aND6EhWw4EfdZCSt0tGQWceB9tXCKVlgc3/TXdTOB9zuyoZxkmQ6uv
rV2ffxf2VLwmR6UJSXsAz2Pd9eWJmnH+QmZPMXhOVFCMRTHTsRfAeyLW+q2xVr/r
c1nV/9PzPP29GSYVb54Fs7of2oHUuBOWp3+2oRljPeoz0SEBG/Q0TdmBqfYTol9r
GapIcROc1qg9oHV6dmQMTAkx3+Io8zlbDp3Xu2+QagtCS+94DcH9Yjh8ggM6hohX
2ofP6HQUw4TLHVTLI0iMc3MJcEZ88voQbHWKT9fYniQjKBESU21IErKT3YWP2OAo
c5RR44gCmE+r14mHCktOLLQrR6sBABEBAAGJAiUEGAECAA8FAldEmzkCGwwFCQPC
ZwAACgkQPVkZtEhFfuDc6g/+PFjc5b156juNTuSyU08/lP58vOWBeL9YxAFcRKVq
mqFQUZSs5xkEA0j0PuueHdPAmLPXO0qE3TRHHcMojI6bpnqPVR6bBHKSE+OUk2yZ
q0LRfsn/Qmn2waIDOKOxA82MAAETiG7/E+obUVGAN2E1fZ30A45eGgqcfr85QYWn
oAWhsenE4gO54ltA8YWyjYvAn+XZ7IegRQE8/Tn2sbKkg8wq/9xnX5rzwVFwanB9
Oki7djcTKeLtt9xBV6l00icB1bYAPiHGXJiBQElDXiUGqyS/yzCL0El6mbNY0NlO
kSWtfu2Rtg5pofJqs4MA8bF84//9JCHRDpQ35uBtfjhlDP0oCA1uKGM85MbBzpv/
Be12Fewqy+ofTLZZk1xPUWIUmRggw0cV2PrJpCk0Dg0E4odRI0niY4OjuA8+yiSS
AJAgYs/fSKLqbF7XMLo9yJ56Z/jrDLqO3mYUhcRAg7EiopE/VdEbFppIIu5rZ1hO
PgzFIwsuI+wrcVI1oYYbVD3gIUdJhe8QSWJ7+x1o5vQRJk0vNwh6eBJwCCeQOKn6
ikXKkkQbc3vsff7wLByyyIKs6zR8aVJk3AHs1QNBusqGZTmBRCrNHg9OJOWfrXSL
WL5zxHbfX+wspgLDLIPYFKlhSEqnfbhrK/17GAd/YF9O7TFy6FzLprmJgp9TaGao
iSI=
=i9Ui
-----END PGP PUBLIC KEY BLOCK-----
EOF

if [ "${RELEASE}" == "trusty" ]; then
    if [ ! -f /etc/apt/sources.list.d/openjdk-r-ppa-trusty.list ]
    then
        sudo apt-add-repository -y ppa:openjdk-r/ppa
    fi
fi

apt-get -yq update

mkdir -p /share

ln -sf /share/reference /reference

apt-get -yq install openjdk-8-{jre,jdk}{,-headless} bazel make libbz2-dev libcurl3-dev libxml2-dev wamerican python-requests git-flow autoconf automake libtool google-perftools libgoogle-perftools-dev

if [ "${RELEASE}" == "xenial" ]; then
	apt-get -yq install python-boto3 clang libomp-dev
fi

apt-get -yq purge nginx nginx-common nginx-core

curl https://raw.githubusercontent.com/bobthecow/git-flow-completion/master/git-flow-completion.bash > /etc/bash_completion.d/git-flow

# ssh agent forwarding support
(
cat <<-"EOM"

# ssh agent forwarding support
cat <<-EOF > ~/.ssh/current-agent
  export SSH_AUTH_SOCK='$SSH_AUTH_SOCK'
  export SSH_CLIENT='$SSH_CLIENT'
  export SSH_CONNECTION='$SSH_CONNECTION'
  export SSH_TTY='$SSH_TTY'
EOF
EOM
) >> ~/.bash_profile

(
cat <<-"EOM"

# refresh the ssh agent
export PROMPT_COMMAND='. ~/.ssh/current-agent'
EOM
) >> ~/.bashrc


# Platform specific stuff

if [ ${PLATFORM} == "AWS" ]; then
    echo 'fs-831e9d2a.efs.us-west-2.amazonaws.com:/       /share  nfs4    nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2,_netdev,nofail        0       0' >> /etc/fstab

    DEST='/usr/local/sbin'
    aws s3 cp s3://sg-other/spiral-repo/ec2mgr/ec2mgr ${DEST}/ec2mgr
    chmod a+x ${DEST}/ec2mgr
    aws s3 cp s3://sg-other/spiral-repo/ec2mgr/ec2mgr-warn.sh /etc/profile.d/
    echo "@reboot      root ${DEST}/ec2mgr boot" > /etc/cron.d/ec2mgr
    echo "0,30 * * * *   root ${DEST}/ec2mgr poll" >> /etc/cron.d/ec2mgr

else
    echo "Preparing SSD..."
    tune2fs -m 0 /dev/sdb1

    mkdir -p /share

    # Azure
    echo
    echo "********************************************************************************"
    echo "You will need to mount the /share directory manually. Run the following as root,"
    echo "substituting the real share password:"
    echo
    echo "    echo '//sgdevshare.file.core.windows.net/sgdevshare /share cifs vers=3.0,username=sgdevshare,password=ACTUAL_PASSWORD_HERE,dir_mode=0777,file_mode=0777,sec=ntlmssp 0 0' >> /etc/fstab"
    echo
    echo "    mount /share"
    echo
    echo "********************************************************************************"
fi
