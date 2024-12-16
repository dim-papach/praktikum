# Table of Contents

1. [Quick installation](#org3da06c5)
2. [Quick Guides](#org98e429a)

This project is based on the work of (González-Gaitán, S and {de~Souza}, R S and Krone-Martins, A and Cameron, E and Coelho, P and Galbany, L and Ishida, E E O and {COIN collaboration}, 2019) and (Smole, Majda and Rino-Silvestre, João and González-Gaitán, Santiago and Stalevski, Marko, 2023).

It uses the R-INLA package to fill missing values in FITS images of astrophysical objects. Our main focus is to test the accuracy and efficiency of this method.

To make sure you don&rsquo;t have a dependency problem, follow the instructions

# Quick installation

- Make sure you have `Nix` installed

- Clone the repository
  
  ```shell
  git clone # clone the repository
  
  cd # move into the directory
  ```

- Build the Nix environment
  
  ```shell
  nix-build
  ```

- Initiate and use the environment
  
  ```shell
  nix-shell
  ```
  
  This opens a bash shell in you current directory and you can run your scripts
  
  ---

# Quick Guides

- [Git](./guides/git.md)
- [Nix](./guides/nix.md)
