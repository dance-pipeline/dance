if __name__ == "__main__":
    with open("part.in", "w") as dict_file:
        methods = ["xtb", "am1"]
        for m in methods:
            dict_file.write(f"{m},{m}.sdf,Energy {m}\n")

