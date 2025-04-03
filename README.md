# Airline-Graphs

## Build Instructions
### 1. **Clone or Download the Source**

**Note**: The easiest way to build is through github codespaces

Make sure you have the source code in your project directory(github codespaces handles this for you).

It can be cloned using
```bash
git clone https://github.com/wba6/Airline-Graphs.git
```
Move into the project directory
```bash
cd ./Airline-Graphs
```

### 3. **Create a Build Directory**
It’s a good practice to create a separate directory for out-of-source builds. From your project root directory, run:
```bash
mkdir build
cd build
```

### 4. **Configure the Build**
Run CMake to configure the project. You can do this with:
```bash
cmake ..
```
This command tells CMake to generate the build files in the current directory (`build`) using the `CMakeLists.txt` in the parent directory.

### 5. **Build the Project**
Once configuration is complete, compile the project by running:
```bash
cmake --build .
```
Alternatively, if you prefer using `make`:
```bash
make
```
This will compile the project and produce an executable named `Airline-Paths`.

### 6. **Run the Executable**
After the build completes, you can run the executable. For example, from the `build` directory:
```bash
./Airline-Paths
```

**Note:** All files are relative to the executable file

### Troubleshooting
- **Compiler Issues:** Ensure that you have a C++ compiler that supports C++17 and that it is properly set up in your environment.
