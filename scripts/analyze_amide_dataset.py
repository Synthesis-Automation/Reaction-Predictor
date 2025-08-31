#!/usr/bin/env python3
"""
Analyze the amide formation dataset to extract analytics for recommendation system.
"""

import json
import os
from collections import Counter

def analyze_amide_dataset():
    """Analyze amide formation dataset and create analytics summary."""
    
    # Counters for different components
    reagent_counts = Counter()
    solvent_counts = Counter()
    base_counts = Counter()
    catalyst_counts = Counter()
    temp_values = []
    time_values = []
    yield_values = []
    
    dataset_path = '../data/reaction_dataset/amide-formation-2021-2024.jsonl'
    if not os.path.exists(dataset_path):
        dataset_path = 'data/reaction_dataset/amide-formation-2021-2024.jsonl'
    
    if not os.path.exists(dataset_path):
        print(f"Dataset not found at {dataset_path}")
        return
    
    print(f"Analyzing dataset: {dataset_path}")
    
    with open(dataset_path, 'r', encoding='utf-8') as f:
        reaction_count = 0
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            try:
                reaction = json.loads(line)
                reaction_count += 1
                
                # Analyze reagents (coupling agents)
                for reagent in reaction.get('reagents', []):
                    name = reagent.get('name', 'Unknown')
                    role = reagent.get('role', '')
                    
                    reagent_counts[name] += 1
                    
                    # Separate bases
                    if role == 'BASE':
                        base_counts[name] += 1
                
                # Analyze solvents
                for solvent in reaction.get('solvents', []):
                    name = solvent.get('name', 'Unknown')
                    solvent_counts[name] += 1
                
                # Analyze catalysts
                catalyst_info = reaction.get('catalyst', {})
                for catalyst in catalyst_info.get('full_system', []):
                    name = catalyst.get('name', 'Unknown')
                    catalyst_counts[name] += 1
                
                # Collect conditions
                conditions = reaction.get('conditions', {})
                temp = conditions.get('temperature_c')
                if temp is not None:
                    temp_values.append(temp)
                    
                time = conditions.get('time_h')
                if time is not None:
                    time_values.append(time)
                    
                yield_val = conditions.get('yield_pct')
                if yield_val is not None:
                    yield_values.append(yield_val)
                    
            except json.JSONDecodeError:
                continue
    
    print(f"\nAnalyzed {reaction_count} reactions")
    
    # Create analytics summary
    analytics = {
        "reaction_type": "Amide Formation",
        "dataset_info": {
            "total_reactions": reaction_count,
            "source_file": "amide-formation-2021-2024.jsonl"
        },
        "top": {
            "reagents": [
                {"name": name, "count": count, "percentage": count/reaction_count*100}
                for name, count in reagent_counts.most_common(20)
            ],
            "solvents": [
                {"name": name, "count": count, "percentage": count/reaction_count*100}
                for name, count in solvent_counts.most_common(15)
            ],
            "bases": [
                {"name": name, "count": count, "percentage": count/reaction_count*100}
                for name, count in base_counts.most_common(10)
            ],
            "catalysts": [
                {"name": name, "count": count, "percentage": count/reaction_count*100}
                for name, count in catalyst_counts.most_common(10)
            ]
        },
        "conditions": {
            "temperature": {
                "avg": sum(temp_values)/len(temp_values) if temp_values else None,
                "min": min(temp_values) if temp_values else None,
                "max": max(temp_values) if temp_values else None,
                "samples": len(temp_values)
            },
            "time": {
                "avg": sum(time_values)/len(time_values) if time_values else None,
                "min": min(time_values) if time_values else None,
                "max": max(time_values) if time_values else None,
                "samples": len(time_values)
            },
            "yield": {
                "avg": sum(yield_values)/len(yield_values) if yield_values else None,
                "min": min(yield_values) if yield_values else None,
                "max": max(yield_values) if yield_values else None,
                "samples": len(yield_values)
            }
        }
    }
    
    # Print summary
    print("\n=== TOP REAGENTS (Coupling Agents) ===")
    for item in analytics["top"]["reagents"][:10]:
        print(f"  {item['name']}: {item['count']} times ({item['percentage']:.1f}%)")
    
    print("\n=== TOP SOLVENTS ===")
    for item in analytics["top"]["solvents"][:10]:
        print(f"  {item['name']}: {item['count']} times ({item['percentage']:.1f}%)")
    
    print("\n=== TOP BASES ===")
    for item in analytics["top"]["bases"][:5]:
        print(f"  {item['name']}: {item['count']} times ({item['percentage']:.1f}%)")
    
    print("\n=== TOP CATALYSTS ===")
    for item in analytics["top"]["catalysts"][:5]:
        print(f"  {item['name']}: {item['count']} times ({item['percentage']:.1f}%)")
    
    # Save analytics to file
    analytics_dir = '../data/analytics/AmideFormation'
    if not os.path.exists(analytics_dir):
        analytics_dir = 'data/analytics/AmideFormation'
    
    os.makedirs(analytics_dir, exist_ok=True)
    
    output_file = os.path.join(analytics_dir, 'latest.json')
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(analytics, f, indent=2, ensure_ascii=False)
    
    print(f"\nAnalytics saved to: {output_file}")
    return analytics

if __name__ == "__main__":
    analyze_amide_dataset()
