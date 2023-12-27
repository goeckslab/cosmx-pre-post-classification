from enum import Enum


class Treatment(Enum):
    ON_PROGRESSION = "ON"
    PRE_TREATMENT = "PRE"


BIOPSY_MAPPINGS = {
    "v1": {
        "265303": Treatment.ON_PROGRESSION.value,
        "321955": Treatment.ON_PROGRESSION.value,
        "322078": Treatment.ON_PROGRESSION.value,
        "272830": Treatment.PRE_TREATMENT.value,
        "303148": Treatment.ON_PROGRESSION.value,
        "272840": Treatment.PRE_TREATMENT.value,
        "321920": Treatment.ON_PROGRESSION.value,
    },
    "v2": {
        "9_2_1": Treatment.PRE_TREATMENT.value,
        "9_2_2": Treatment.ON_PROGRESSION.value,
        "9_3_1": Treatment.PRE_TREATMENT.value,
        "9_3_2": Treatment.ON_PROGRESSION.value,
        "9_14_1": Treatment.PRE_TREATMENT.value,
        "9_14_2": Treatment.ON_PROGRESSION.value,
        "9_15_1": Treatment.PRE_TREATMENT.value,
        "9_15_2": Treatment.ON_PROGRESSION.value,
    },
    "v3": {
        "209184_1": Treatment.PRE_TREATMENT.value,
        "209184_2": Treatment.ON_PROGRESSION.value,
        "312158_1": Treatment.PRE_TREATMENT.value,
        "312158_2": Treatment.ON_PROGRESSION.value,
        "475139_1": Treatment.PRE_TREATMENT.value,
        "475139_2": Treatment.ON_PROGRESSION.value,
        "484265_1": Treatment.PRE_TREATMENT.value,
        "484265_2": Treatment.ON_PROGRESSION.value,
        "566112_1": Treatment.PRE_TREATMENT.value,
        "566112_2": Treatment.ON_PROGRESSION.value,
        "657907_1": Treatment.PRE_TREATMENT.value,
        "657907_2": Treatment.ON_PROGRESSION.value,
        "672431_1": Treatment.PRE_TREATMENT.value,
        "672431_2": Treatment.ON_PROGRESSION.value,
        "798747_1": Treatment.PRE_TREATMENT.value,
        "798747_2": Treatment.ON_PROGRESSION.value,
        "894222_1": Treatment.PRE_TREATMENT.value,
        "894222_2": Treatment.ON_PROGRESSION.value,
    }
}

