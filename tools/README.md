# ProtFlow Tools Directory

This folder is reserved for optional external local bio-computation tools.

- **P2Rank**: Extracted to `tools/p2rank_2.5.1/` or `tools/p2rank/`
- **Java JRE**: Extracted to `tools/jre/` (or via system `JAVA_HOME`)
- **MMseqs2**: Extracted to `tools/mmseqs/`

ProtFlow automatically detects these if present; if not installed, ProtFlow safely falls back to built-in geometric algorithms and remote APIs.
Heavy binaries and runtime caches in this directory are excluded from Git by `.gitignore`.
