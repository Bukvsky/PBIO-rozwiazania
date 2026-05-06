import random
import os


def generate_sequence(length: int) -> str:
    nukleotydy = ["A", "C", "G", "T"]
    seq = ""
    for _ in range(length):
        seq += random.choice(nukleotydy)
    return seq


def calculate_stats(sequence: str) -> dict:
    n = len(sequence)
    if n == 0:
        return {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0, "GC": 0.0}

    a = sequence.count("A")
    c = sequence.count("C")
    g = sequence.count("G")
    t = sequence.count("T")

    stats = {
        "A": round(a / n * 100, 2),
        "C": round(c / n * 100, 2),
        "G": round(g / n * 100, 2),
        "T": round(t / n * 100, 2),
        "GC": round((g + c) / n * 100, 2),
    }
    return stats


def insert_name(sequence: str, name: str) -> str:
    pos = random.randint(0, len(sequence))
    return sequence[:pos] + name.lower() + sequence[pos:]


def format_fasta(seq_id: str, description: str,
                 sequence: str, line_width: int = 80) -> str:
    if description:
        header = ">" + seq_id + " " + description
    else:
        header = ">" + seq_id

    lines = []
    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i + line_width])

    return header + "\n" + "\n".join(lines) + "\n"


def validate_positive_int(prompt: str,
                          min_val: int = 1,
                          max_val: int = 100_000) -> int:
    while True:
        wartosc = input(prompt)
        try:
            n = int(wartosc)
            if n < min_val or n > max_val:
                print(f"Błąd: wartość musi być liczbą całkowitą z zakresu [{min_val}, {max_val}].")
                continue
            return n
        except ValueError:
            print(f"Błąd: wartość musi być liczbą całkowitą z zakresu [{min_val}, {max_val}].")


def get_id():
    while True:
        seq_id = input("Podaj ID sekwencji: ").strip()
        if seq_id == "":
            print("Błąd: ID nie może być puste.")
            continue
        if any(c.isspace() for c in seq_id):
            print("Błąd: ID nie może zawierać białych znaków.")
            continue
        return seq_id


def find_motif(sequence: str, motif: str) -> list:
    motif = motif.upper()
    pozycje = []
    for i in range(len(sequence) - len(motif) + 1):
        if sequence[i:i + len(motif)] == motif:
            pozycje.append(i + 1)
    return pozycje


def complement(sequence: str) -> str:
    pary = {"A": "T", "T": "A", "C": "G", "G": "C"}
    out = ""
    for n in sequence:
        out += pary[n]
    return out


def reverse_complement(sequence: str) -> str:
    return complement(sequence)[::-1]


def transcribe(sequence: str) -> str:
    return sequence.replace("T", "U")


def find_orfs(sequence: str, min_length: int) -> list:
    stop_kodony = ["TAA", "TAG", "TGA"]
    orfy = []

    for ramka in range(3):
        i = ramka
        while i < len(sequence) - 2:
            kodon = sequence[i:i + 3]
            if kodon == "ATG":

                j = i + 3
                while j < len(sequence) - 2:
                    if sequence[j:j + 3] in stop_kodony:
                        dlugosc = j + 3 - i
                        if dlugosc >= min_length:
                            orfy.append((i + 1, j + 3, dlugosc))
                        i = j + 3
                        break
                    j += 3
                else:

                    break
            else:
                i += 3
    return orfy


def main():
    print(" Generator sekwencji DNA -> FASTA ")

    length = validate_positive_int("Podaj długość sekwencji: ")
    seq_id = get_id()
    description = input("Podaj opis sekwencji: ").strip()
    name = input("Podaj imię: ").strip()
    while name == "" or not name.isalpha():
        print("Błąd: imię musi zawierać tylko litery.")
        name = input("Podaj imię: ").strip()

    sequence = generate_sequence(length)
    stats = calculate_stats(sequence)

    seq_with_name = insert_name(sequence, name)

    filename = f"{seq_id}.fasta"
    fasta_record = format_fasta(seq_id, description, seq_with_name)

    with open(filename, "w") as f:
        f.write(fasta_record)

    print(f"\nSekwencja zapisana do pliku: {filename}")
    print(f"\nStatystyki sekwencji (n={length}):")
    print(f"  A: {stats['A']:.2f}%")
    print(f"  C: {stats['C']:.2f}%")
    print(f"  G: {stats['G']:.2f}%")
    print(f"  T: {stats['T']:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")


    # 1. wyszukiwanie motywow
    print("\n Wyszukiwanie motywu ")
    motif = input("Podaj motyw do wyszukania (lub Enter aby pominąć): ").strip().upper()
    if motif != "":
        # walidacja motywu
        if all(c in "ACGT" for c in motif):
            wyniki = find_motif(sequence, motif)
            if wyniki:
                print(f"Motyw '{motif}' znaleziony na pozycjach: {wyniki}")
                print(f"Liczba wystąpień: {len(wyniki)}")
            else:
                print(f"Motyw '{motif}' nie został znaleziony.")
        else:
            print("Motyw zawiera niedozwolone znaki - pominięto.")

    # 2. komplementarna i odwrotnie komplementarna 
    comp = complement(sequence)
    rev_comp = reverse_complement(sequence)
    with open(filename, "a") as f:
        f.write(format_fasta(seq_id + "_complement",
                             "nic komplementarna", comp))
        f.write(format_fasta(seq_id + "_revcomp",
                             "nic odwrotnie komplementarna", rev_comp))
    print("\nDopisano nici komplementarne do pliku FASTA.")

    # 3. transkrypcja
    mrna = transcribe(sequence)
    with open(filename, "a") as f:
        f.write(format_fasta(seq_id + "_mRNA",
                             "transkrypt mRNA", mrna))
    print("Dopisano transkrypt mRNA do pliku FASTA.")

  
    print("\n Identyfikacja ORF ")
    min_len = validate_positive_int("Podaj minimalną długość ORF (nt): ", 3, 100_000)
    orfy = find_orfs(sequence, min_len)
    if orfy:
        print(f"Znaleziono {len(orfy)} ORF:")
        for start, stop, dl in orfy:
            print(f"  pozycja {start}-{stop}, długość: {dl} nt")
    else:
        print("Nie znaleziono żadnych ORF spełniających kryterium.")


if __name__ == "__main__":
    main()