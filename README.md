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

If you run a Linux distribution the easiest way to run create and run the containers is distrobox, which takes care of linking your physical resources with the virtual ones of the container and mapping your home directory correctly. Distrobox should be already available to install in the package manager of your distribution, if not already installed by default (i.e. Fedora). If this is not the case, please follow the official documentation: [https://distrobox.it/#installation](https://distrobox.it/#installation).

After installing distrobox create the desider container with:

```
$ distrobox create -n groot -i docker.io/alexo91/groot:XXXX
```
Replace `XXXX` with either `ubuntu` or `almalinux` tag to choose your preferred version. If in doubt we suggest you to use the ubuntu version.
If you need root permission add `--root` at the end of the command.

Then to enter the container simply run:

```
$ distrobox enter groot
```
Once inside you should run

```
$ source /opt/geant4/bin/geant4.sh
$ source /opt/root/bin/thisroot.sh
$ export PATH="/opt/GROOT/bin:$PATH"
```
Then you can run GROOT, GEANT4 and ROOT from everywhere. 

Other than distrobox, there are other utilities available to achieve the same end, such as [x11docker](https://github.com/mviereck/x11docker). If you want to use one of those, please follow the respective documentation.

## macOS

If you don't have much experience with docker, we suggest you to install Docker Desktop to manage and run your containers, by following the simple guide here [https://docs.docker.com/desktop/install/mac-install/](https://docs.docker.com/desktop/install/mac-install/).

Once docker is installed in your system and you can run the `docker` command from the terminal, you should install XQuartz by either downloading the installer from the [official website](https://www.xquartz.org/) or by using homebrew with `brew install --cask xquartz`.

After installing XQuartz open it and go to Security Settings and ensure that "Allow connections from network clients" is on, then reboot your Mac and start XQuartz again.

You should now allow X11 forwarding for local containers via xhost (**you should run this each time you restart XQuartz!**)

```
$ xhost +localhost
```

Pull the docker container from the registry

```
$ docker pull alexo91/groot:XXXX
```

Replace `XXXX' with either `ubuntu` or `almalinux` tag to choose your preferred version. If in doubt we suggest you to use the ubuntu version.
And then run the container with:

```
$ docker run -rm -e DISPLAY=host.docker.internal:0 -v /tmp/.X11-unix:/tmp/.X11-unix alexo91/groot:XXXX
```
If you need to share a local folder of your mac with the GROOT container please specify it in the command above with the syntax

```
-v /mac/local/path:/path/in/the/container:rw
```
The `rw` at the end stands for "read&write" permissions, for read-only use `ro`.

If you want to open the bash terminal instead of directly running GROOT you can use the following command:

```
$ docker run -rm -it -e DISPLAY=host.docker.internal:0 -v /tmp/.X11-unix:/tmp/.X11-unix alexo91/groot:XXXX bash
```

## Windows

If you don't have much experience with docker, we suggest you to install Docker Desktop to manage and run your containers, by following the guide here [https://docs.docker.com/desktop/install/windows-install/](https://docs.docker.com/desktop/install/windows-install/).

Once docker is installed you need to nstall an X Server (such as Xming, VcXsrv or another X server of choice) and configure it to allow connections from other hosts. For Xming, you can use the “-ac” option to disable access control. From the command prompt run:

```
Xming :0 -ac -multiwindow
```


Pull the docker container from the registry

```
$ docker pull alexo91/groot:XXXX
```

Replace `XXXX' with either `ubuntu` or `almalinux` tag to choose your preferred version. If in doubt we suggest you to use the ubuntu version.
And then run the container with:

```
$ docker run -rm -e DISPLAY=host.docker.internal:0 alexo91/groot:XXXX
```
If you need to share a local folder of your mac with the GROOT container please specify it in the command above with the syntax

```
-v /pc/local/path:/path/in/the/container:rw
```
The `rw` at the end stands for "read&write" permissions, for read-only use `ro`.

If you want to open the bash terminal instead of directly running GROOT you can use the following command:

```
$ docker run -rm -it -e DISPLAY=host.docker.internal:0 alexo91/groot:XXXX bash
```

