# FEDDLib Build System with TriBITS

This directory contains the build system configuration for FEDDLib, which uses **TriBITS** (Tribal Build, Integrate, and Test System) as its CMake-based build framework.
It was copied directly from the TriBITS github (see below) at commit 02fb9c95716b1a13ba7aa841e82cbe9ac93a2c34 and trimmed down to suit the needs of the FEDDLib.

FEDDLib is organized into the following TriBITS packages:

```
FEDDLib/
├── core/           # Core finite element and domain decomposition functionality
├── problems/       # Problem-specific implementations and solvers
└── amr/           # Adaptive mesh refinement capabilities
```

### Package Dependencies

- **core**: Base package with no dependencies
- **problems**: Depends on `core`
- **amr**: Depends on `core` and `problems`

## Key Configuration Structure

- `../CMakeLists.txt` - Main project configuration
- `../ProjectName.cmake` - Project name and basic settings
- `../PackagesList.cmake` - Package definitions and classifications
- `../TPLsList.cmake` - Third-party library definitions
- `../Version.cmake` - Version information
- `tribits/` - TriBITS framework (subset of full TriBITS)
- `TPLs/FindTPLTrilinos.cmake` - Trilinos integration
- `TPLs/FindTPLAceGENInterface.cmake` - AceGEN interface (optional)

## TriBITS Resources

- **Main Website**: [https://tribits.org](https://tribits.org)
- **GitHub Repository**: [https://github.com/TriBITSPub/TriBITS](https://github.com/TriBITSPub/TriBITS)
- **Example Projects**: [https://github.com/TriBITSPub/TriBITS/tree/master/tribits/examples](https://github.com/TriBITSPub/TriBITS/tree/master/tribits/examples)
