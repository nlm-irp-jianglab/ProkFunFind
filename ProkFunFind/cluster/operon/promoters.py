from typing import Any, Dict, List, Union, Set, Tuple  # pylint: disable=unused-import

class Promoter:
    """ Contains all the relevant info and helpers for promoters """

    def __init__(self, gene_name: str, start: int, end: int, seq: str = "") -> None:
        self.gene_name = str(gene_name)
        self.start = int(start)
        self.end = int(end)
        self.seq = seq

    def get_id(self) -> str:
        """ Returns the id of the promoter """
        return self.gene_name

    def get_gene_names(self) -> List[str]:
        """ Returns a list of gene names attached to this promoter """
        return [self.gene_name]

    def __len__(self) -> int:
        if self.seq is None:
            raise ValueError(
                "Requesting length of a promoter sequence which hasn't been set")
        return len(self.seq)

    def __str__(self) -> str:
        return "Promoter(%r, %d, %d)" % (self.get_id(), self.start, self.end)

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other: Any) -> bool:
        return (isinstance(other, Promoter)
                and self.gene_name == other.gene_name
                and self.start == other.start
                and self.end == other.end
                and self.seq == other.seq)

    def to_json(self) -> Dict[str, Any]:
        """ Serialises a Promoter to a dictionary for use in JSON """
        return {"gene": self.gene_name, "start": self.start, "end": self.end,
                "seq": str(self.seq), "type": "Promoter"}

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "Promoter":
        """ Deserialises a Promoter from a dictionary """
        return Promoter(str(json["gene"]), int(json["start"]), int(json["end"]), seq=str(json["seq"]))
