from enum import Enum


class Treatment(Enum):
    ON_PROGRESSION = "ON"
    PRE_TREATMENT = "PRE"


BIOPSY_PRE_OR_POST = {
    "265303": Treatment.ON_PROGRESSION.value,
    "321955": Treatment.ON_PROGRESSION.value,
    "322078": Treatment.ON_PROGRESSION.value,
    # "265607": Treatment.PRE_TREATMENT.value,
    # "322118": Treatment.ON_PROGRESSION.value,
    "272830": Treatment.PRE_TREATMENT.value,
    "303148": Treatment.ON_PROGRESSION.value,
    "272840": Treatment.PRE_TREATMENT.value,
    "321920": Treatment.ON_PROGRESSION.value,
    # "322091": Treatment.ON_PROGRESSION.value
}
