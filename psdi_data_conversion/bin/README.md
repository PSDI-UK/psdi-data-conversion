# Converter Binaries

This folder contains precompiled binaries for converters, targeting the most common platforms. These binaries are intended to be runnable without installing any additional dependencies such as LAPACK.

## Structure

The files here are structured as: `<os>/<plugin>`, where `<os>` is a folder named either "linux" or "mac" and `<plugin>` is a plugin binary - currently either "atomsk" or "c2x". The plugin binaries do not have to be named matching the plugin, and plugins aren't limited to a single binary; the use of these binaries is determined by their execution scripts (see the `psdi_data_conversion/scripts` folder).

This folder also includes licenses for the Atomsk and c2x converters. These are both licensed under the GPLv3, which allows us to provide precompiled binaries here provided that we include a copy of their license, and this serves to meet that obligation. Any future plugins added with a similar requirement should have their licenses included here as well.

## Obtaining binaries

### Atomsk

We currently package version 0.13.1 of Atomsk. When this is updated, the version here should be updated as well.

Atomsk provides binary distributions from its website at https://atomsk.univ-lille.fr/, on the [Downloads page](https://atomsk.univ-lille.fr/dl.php). For the Linux binary, we use the "Linux amd64" package, which includes the binary as `atomsk` within the self-named folder at the root of the tarball. This file can be copied to `linux/atomsk` to update it, replacing the existing binary. It may be necessary to manually enable execution permission (`chmod +x atomsk`) depending on how the file is copied. This binary is statically-linked and should be portable across many Linux distributions.

Atomsk does not provide a precompiled MacOS binary, so we set up a fork of the Atomsk source repository at https://github.com/PSDI-UK/atomsk with a GitHub workflow to compile it at https://github.com/PSDI-UK/atomsk/actions/workflows/makefile.yml. This workflow uses GitHub's `macos-latest` runners (currently macOS 26.6.2) to compile the binary, which is uploaded as an artifact. This artifact can be downloaded, and the `atomsk` binary copied to `macos/atomsk` with execution permission enabled (`chmod +x atomsk`) to update it. This binary is dynamically linked to the BLAS and LAPACK packages which come built-in on this OS, and so may not be compatible with different OS versions.

We also investigated compiling our own Linux binary through the same repo, in case the Atomsk site becomes unavailable in the future. Unfortunately, we were unable to recreate their method of generating a statically-linked binary (using their provided Makefile for a static binary results in a segmentation fault on execution of the generated binary). This code is left in the repository in case it is needed in the future.

### c2X

c2x provides binary distributions from its website at https://www.c2x.org.uk/, on the [Downloads page](https://www.c2x.org.uk/downloads/). We use the Linux distribution and the "MacOS, Intel, 2.42a, with spglib 2.0.2" (or similar in the future) distribution, which may be downloaded and then copied to `<os>/c2x` with execution permission enabled (`chmod +x c2x`) to update them. Both of these binaries appear to work well across different OS versions.

In case the c2x site becomes unavailable in the future or we need to recompile the binaries for future OS versions, we can pursue a similar solution as with Atomsk, setting up a fork of the repository with a workflow to build binaries for the desired distributions.
