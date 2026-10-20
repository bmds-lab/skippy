from pathlib import Path

project = "Skippy"
extensions = ["breathe"]

breathe_projects = {
    "Skippy": str(Path(__file__).parent / "_build" / "xml")
}
breathe_default_project = "Skippy"
