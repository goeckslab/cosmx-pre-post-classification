from enum import Enum


class Treatment(Enum):
    ON_PROGRESSION = "ON"
    PRE_TREATMENT = "PRE"


BIOPSY_MAPPINGS = {
    "265303": Treatment.ON_PROGRESSION.value,
    "321955": Treatment.ON_PROGRESSION.value,
    "322078": Treatment.ON_PROGRESSION.value,
    "272830": Treatment.PRE_TREATMENT.value,
    "303148": Treatment.ON_PROGRESSION.value,
    "272840": Treatment.PRE_TREATMENT.value,
    "321920": Treatment.ON_PROGRESSION.value,
    # "322091": Treatment.ON_PROGRESSION.value

    "9_2_1": Treatment.PRE_TREATMENT.value,
    "9_2_2": Treatment.ON_PROGRESSION.value,
    "9_3_1": Treatment.PRE_TREATMENT.value,
    "9_3_2": Treatment.ON_PROGRESSION.value,
    "9_14_1": Treatment.PRE_TREATMENT.value,
    "9_14_2": Treatment.ON_PROGRESSION.value,
    "9_15_1": Treatment.PRE_TREATMENT.value,
    "9_15_2": Treatment.ON_PROGRESSION.value,

}
