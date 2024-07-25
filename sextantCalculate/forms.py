# forms.py
from django import forms
from .utils import load_star_data

# Load star data and create choices with star names
star_data = load_star_data('sextantCalculate/NavigationalStars.csv')
STAR_CHOICES = [(i, f"{star[1]} ({i})") for i, star in enumerate(star_data)]

class ObservationForm(forms.Form):
    star_number = forms.ChoiceField(choices=STAR_CHOICES)
    elevation = forms.FloatField()
