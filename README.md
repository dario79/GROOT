# GROOT

GROOT makes use of the validated tracking GEANT4 libraries and the n-body event generator of ROOT in order to provide a fast, realiable and complete simulation tool to be used for nuclear physics experiments.


# Installation

## Precompiled binaries

GROOT is available, alongside the required GEANT4 and ROOT binaries and libraries, as a precompiled binary. To download the archive please refer to the realease section of this repository. At the moment it built for only the latest Ubuntu LTS (24.04) and the latest Alma Linux 9.
Before running the program you must install the required dependencies. Below there are the terminal commands as one-liners. Please be aware that you need to execute them as root or using sudo (with an account in the sudoers group).

- Ubuntu 24.04
  ```
  # apt install qtbase5-dev libqt5opengl5-dev libxmu-dev \
    libgl1-mesa-dev qt3d5-dev libxerces-c-dev qtbase5-dev \
    libssl-dev git libx11-dev \
    libxext-dev libxft-dev libxpm-dev python3 libtbb-dev \
    gfortran libpcre3-dev \
    libglu1-mesa-dev libglew-dev libftgl-dev \
    libfftw3-dev libcfitsio-dev libgraphviz-dev \
    libavahi-compat-libdnssd-dev libldap2-dev \
    python3-dev python3-numpy libxml2-dev libkrb5-dev \
    libgsl-dev qtwebengine5-dev nlohmann-json3-dev libmysqlclient-dev
  ```
- Alma Linux 9
  ```
  # dnf install xerces-c-devel qt5-devel \
    libX11-devel libXpm-devel libXmu-devel libXft-devel libXext-devel python openssl-devel \
    xrootd-client-devel xrootd-libs-devel \
    mesa-libGL-devel mesa-libGLU-devel glew-devel ftgl-devel mysql-devel \
    fftw-devel cfitsio-devel graphviz-devel libuuid-devel \
    avahi-compat-libdns_sd-devel openldap-devel python-devel python3-numpy \
    libxml2-devel gsl-devel readline-devel qt5-qtwebengine-devel \
  ```

After installing these packages, download the archive corresponding to your distribution and open a terminal in your download folder. Then unpack it in your home directory using the following command:

```
$ tar -xvf groot-XXXX.tar.gz -C $HOME
```
Please replace the `XXXX` with the name of your distro.


## Containers (Docker, Podman, LXC, etc.) 

For all the people which use a different distribution or a different OS (Windows, macOS, etc.) two OCI images, either with Ubuntu 24.04 LTS or Alma Linux 9, are available. In this images you will find GROOT installed alongside GEANT4 and ROOT. You can use your preferred container system to run those, but the most famous are probably Docker and Podman. We will explain how to use Docker with each operating system. Please note that GROOT is mainly run as a GUI application, however Docker was not developed to this aim. Therefore we would need to use some workaround to achieve the final result.

### Linux

If you run a Linux distribution the easiest way to run create and run the containers is distrobox.

## macOS

## Windows



