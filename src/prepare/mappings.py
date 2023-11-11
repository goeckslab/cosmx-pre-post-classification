from enum import Enum


class Treatment(Enum):
    POST_TREATMENT = "POST"
    PRE_TREATMENT = "PRE"


BIOPSY_PRE_OR_POST = {
    "265303": Treatment.POST_TREATMENT.value,
    "321955": Treatment.POST_TREATMENT.value,
    "322078": Treatment.POST_TREATMENT.value,
    "265607": Treatment.PRE_TREATMENT.value,
    "322118": Treatment.POST_TREATMENT.value,
    "272830": Treatment.PRE_TREATMENT.value,
    "303148": Treatment.POST_TREATMENT.value,
    "272840": Treatment.POST_TREATMENT.value,
    "321920": Treatment.POST_TREATMENT.value,
    "322091": Treatment.POST_TREATMENT.value
}
