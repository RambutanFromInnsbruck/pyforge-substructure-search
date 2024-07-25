from pydantic import BaseModel


class Molecule(BaseModel):
    molecule_id: int
    molecule_structure: str
