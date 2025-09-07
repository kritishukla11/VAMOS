# Contributing to VAMOS

1. Create a virtual environment and install in editable mode:

   ```bash
   python -m venv .venv
   source .venv/bin/activate  # or .venv\Scripts\activate on Windows
   pip install -e .
   ```

2. Run the CLI:

   ```bash
   vamos map
   vamos cluster
   vamos ml
   vamos stats
   ```

3. Run tests:

   ```bash
   pytest -q
   ```
