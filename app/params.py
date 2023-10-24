github_link = 'https://github.com/antigenomics/tcr-covid-classifier'
zenodo_link = 'https://zenodo.org/record/8362803'
preprint_link = 'https://docs.google.com/document/d/1sMOrdsspdResrqM3fX3lY3Z-Bc6cOXh5dULpZtWT0os/edit#heading=h.vbacfjakxjr'
button_style = 'display: inline-block; padding: 12px 20px; background-color: #ddece4; color: #237c54; text-align: center; text-decoration: none; font-size: 16px; border-radius: 4px;'
abstract_text = '''
* a ML approach applied to two large datasets 
* hundreds of healthy donor/COVID-19 patients repertoires
* inferred TCR biomarkers that were induced by SARS-CoV-2'''
desc_text = f"""
* accomodates the COVID-19 associated clonotypes found in the [study]({preprint_link}).
* HLA info available
* cluster logo can be assessed
* search for similar clonotypes using regex or Levenstein distance 
* visualize the closest clusters. \n
    """
css = '''
<style>
    .stTabs [data-baseweb="tab-list"] button [data-testid="stMarkdownContainer"] p {
    font-size: 1.4rem;
    }
</style>
'''