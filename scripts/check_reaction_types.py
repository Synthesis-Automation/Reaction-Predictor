import json

types = set()
with open('data/reaction_dataset/amide-formation-2021-2024.jsonl', 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if line:
            try:
                data = json.loads(line)
                rtype = data.get('reaction_type', 'None')
                types.add(rtype)
            except:
                continue

print('Reaction types in dataset:')
for t in sorted(types):
    print(f'  - {t}')
