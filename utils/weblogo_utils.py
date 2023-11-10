import weblogo
import os


def create_web_logo(sequences, output_file=None, out_format=None, title=None, units='probability',
                    alphabet=weblogo.seq.unambiguous_protein_alphabet):
    if out_format is None:
        extension = os.path.splitext(output_file)[1] if output_file is not None else ''
        out_format = extension[1:] if extension else 'png'

    seqs = weblogo.seq.SeqList([weblogo.seq.Seq(s, alphabet) for s in sequences], alphabet)
    seqs.alphabet = alphabet
    data = weblogo.LogoData.from_seqs(seqs)
    options = weblogo.LogoOptions()
    if title is not None:
        options.logo_title = title
    options.unit_name = units
    options.show_fineprint = False

    if out_format == 'png':
        options.resolution = 400.0

    format = weblogo.LogoFormat(data, options)

    formatters = {
        'png': weblogo.eps_formatter,
        'svg': weblogo.svg_formatter
    }

    image = formatters[out_format](data, format)
    if output_file is None:
        return image

    with open(output_file, 'wb') as fout:
        fout.write(image)

