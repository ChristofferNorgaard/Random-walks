{
  description = "Fortran development environment with GPU support";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        
        # Common Fortran packages
        fortranPackages = with pkgs; [
          # Compilers
          gfortran13          # GNU Fortran compiler
          gcc13               # GCC for compatibility
          
          # Build tools
          cmake
          gnumake
          gnuplot
          pkg-config
          
          # Libraries
          scalapack           # Parallel linear algebra
          fftw                # Fast Fourier Transform
          hdf5                # HDF5 file format
          netcdf              # NetCDF file format
          netcdffortran       # Fortran bindings for NetCDF
          
          # MPI for parallel computing
          openmpi
          
          # Debugging and profiling
          gdb
          valgrind
          
          # Documentation
          doxygen
          
          # Version control
          git
        ];
        
        # GPU-specific packages (CUDA)
        gpuPackages = with pkgs; [
          cudatoolkit
          # Note: NVIDIA HPC SDK (nvfortran) not available in nixpkgs
          # You may need to install it separately
        ];

      in
      {
        # Development shell
        devShells.default = pkgs.mkShell {
          name = "fortran-dev";
          
          buildInputs = fortranPackages;
          
          shellHook = ''
            echo "🚀 Fortran Development Environment"
            echo "=================================="
            echo "Compiler: $(gfortran --version | head -n1)"
            echo "MPI: $(mpirun --version | head -n1)"
            echo ""
            echo "Available compilers:"
            echo "  - gfortran: GNU Fortran compiler"
            echo "  - mpifort: MPI Fortran wrapper"
            echo ""
            echo "Key libraries:"
            echo "  - OpenBLAS: Optimized BLAS"
            echo "  - LAPACK: Linear algebra"
            echo "  - FFTW: Fast Fourier Transform"
            echo "  - HDF5: Hierarchical data format"
            echo "  - NetCDF: Network Common Data Form"
            echo ""
            echo "Build tools: cmake, make"
            echo "Debugging: gdb, valgrind"
            echo ""
            
            # Set up environment variables
            export FC=gfortran
            export CC=gcc
            export CXX=g++
            
            # Library paths
            export LIBRARY_PATH=${pkgs.lib.makeLibraryPath fortranPackages}
            export LD_LIBRARY_PATH=${pkgs.lib.makeLibraryPath fortranPackages}
            export PKG_CONFIG_PATH=${pkgs.lib.makeSearchPathOutput "dev" "lib/pkgconfig" fortranPackages}
          '';
        };
        
        # GPU development shell (includes CUDA)
        devShells.gpu = pkgs.mkShell {
          name = "fortran-gpu-dev";
          
          buildInputs = fortranPackages ++ gpuPackages;
          
          shellHook = ''
            echo "🚀 Fortran GPU Development Environment"
            echo "======================================"
            echo "Compiler: $(gfortran --version | head -n1)"
            echo "CUDA: ${pkgs.cudatoolkit.version}"
            echo ""
            echo "⚠️  Note: NVIDIA HPC SDK (nvfortran) must be installed separately"
            echo "    Download from: https://developer.nvidia.com/hpc-sdk"
            echo ""
            echo "Available for GPU computing:"
            echo "  - CUDA Toolkit"
            echo "  - OpenACC (requires NVIDIA HPC SDK)"
            echo ""
            
            export FC=gfortran
            export CC=gcc
            export CXX=g++
            
            # Add CUDA to path
            export CUDA_PATH=${pkgs.cudatoolkit}
            export PATH=${pkgs.cudatoolkit}/bin:$PATH
            export LIBRARY_PATH=${pkgs.lib.makeLibraryPath (fortranPackages ++ gpuPackages)}
            export LD_LIBRARY_PATH=${pkgs.lib.makeLibraryPath (fortranPackages ++ gpuPackages)}
          '';
        };
        
        # Minimal shell for quick compilation
        devShells.minimal = pkgs.mkShell {
          name = "fortran-minimal";
          
          buildInputs = with pkgs; [
            gfortran13
            gnumake
            cmake
          ];
          
          shellHook = ''
            echo "Fortran minimal environment (compiler + build tools)"
            export FC=gfortran
          '';
        };
        
        # Example package: compile a Fortran program
        packages.pareto-sampler = pkgs.stdenv.mkDerivation {
          pname = "pareto-sampler";
          version = "0.1.0";
          
          src = ./.;
          
          nativeBuildInputs = [ pkgs.gfortran13 ];
          
          buildPhase = ''
            gfortran -o pareto_sampler pareto_sampler.f90 -O3
          '';
          
          installPhase = ''
            mkdir -p $out/bin
            cp pareto_sampler $out/bin/
          '';
        };
      }
    );
}
