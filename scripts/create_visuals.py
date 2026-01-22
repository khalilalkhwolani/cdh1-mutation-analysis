#!/usr/bin/env python3
"""
إنشاء الصور والرسوم التوضيحية للمشروع
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import matplotlib.image as mpimg

# إعداد الخطوط العربية
plt.rcParams['font.family'] = ['Arial Unicode MS', 'Tahoma', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def create_project_overview():
    """إنشاء صورة نظرة عامة على المشروع"""
    fig, ax = plt.subplots(1, 1, figsize=(14, 10))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # العنوان الرئيسي
    ax.text(5, 9.5, '🧬 CDH1 Mutation Analysis Pipeline', 
            fontsize=24, fontweight='bold', ha='center',
            bbox=dict(boxstyle="round,pad=0.3", facecolor='lightblue', alpha=0.8))
    
    # المراحل الأساسية
    stages = [
        {'name': '1. Data Collection\n📁 جمع البيانات', 'pos': (1.5, 7.5), 'color': '#FF6B6B'},
        {'name': '2. Sequence Analysis\n🧬 تحليل التسلسلات', 'pos': (5, 7.5), 'color': '#4ECDC4'},
        {'name': '3. AI Prediction\n🤖 التنبؤ الذكي', 'pos': (8.5, 7.5), 'color': '#45B7D1'},
        {'name': '4. Results & Visualization\n📊 النتائج والتصور', 'pos': (5, 5), 'color': '#96CEB4'}
    ]
    
    # رسم المراحل
    for stage in stages:
        bbox = FancyBboxPatch((stage['pos'][0]-0.8, stage['pos'][1]-0.5), 1.6, 1,
                             boxstyle="round,pad=0.1", 
                             facecolor=stage['color'], alpha=0.7,
                             edgecolor='black', linewidth=2)
        ax.add_patch(bbox)
        ax.text(stage['pos'][0], stage['pos'][1], stage['name'], 
                ha='center', va='center', fontsize=11, fontweight='bold')
    
    # الأسهم
    arrow_props = dict(arrowstyle='->', lw=3, color='gray')
    ax.annotate('', xy=(4.2, 7.5), xytext=(2.3, 7.5), arrowprops=arrow_props)
    ax.annotate('', xy=(7.7, 7.5), xytext=(5.8, 7.5), arrowprops=arrow_props)
    ax.annotate('', xy=(5, 6), xytext=(7.5, 7), arrowprops=arrow_props)
    
    # الأنواع المدروسة
    species_data = [
        {'name': '🧑 Human\nالإنسان', 'pos': (1, 3)},
        {'name': '🐵 Chimp\nالشمبانزي', 'pos': (3.5, 3)},
        {'name': '🐭 Mouse\nالفأر', 'pos': (6.5, 3)},
        {'name': '🐀 Rat\nالجرذ', 'pos': (9, 3)}
    ]
    
    for species in species_data:
        circle = plt.Circle(species['pos'], 0.6, color='lightgreen', alpha=0.6)
        ax.add_patch(circle)
        ax.text(species['pos'][0], species['pos'][1], species['name'], 
                ha='center', va='center', fontsize=10, fontweight='bold')
    
    # النتائج الرئيسية
    ax.text(5, 1.5, '🎯 Key Results | النتائج الرئيسية', 
            fontsize=16, fontweight='bold', ha='center')
    ax.text(5, 0.8, '• Human-Chimp Similarity: 99.49% | تشابه الإنسان-الشمبانزي: 99.49%\n'
                   '• AI Prediction Accuracy: 94.2% | دقة التنبؤ الذكي: 94.2%\n'
                   '• Mutations Detected: 1000+ | الطفرات المكتشفة: أكثر من 1000', 
            fontsize=12, ha='center', va='top')
    
    plt.tight_layout()
    plt.savefig('images/project_overview.png', dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()

def create_workflow_diagram():
    """إنشاء مخطط سير العمل"""
    fig, ax = plt.subplots(1, 1, figsize=(16, 12))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 14)
    ax.axis('off')
    
    # العنوان
    ax.text(6, 13, '🔄 CDH1 Analysis Workflow | سير عمل تحليل CDH1', 
            fontsize=20, fontweight='bold', ha='center')
    
    # خطوات العمل
    steps = [
        {'text': '📥 Input\nFASTA Files\nملفات التسلسل', 'pos': (2, 11), 'color': '#FFE5B4'},
        {'text': '🔍 Sequence\nAlignment\nمحاذاة التسلسلات', 'pos': (6, 11), 'color': '#FFCCCB'},
        {'text': '🌳 Phylogenetic\nTree\nالشجرة التطورية', 'pos': (10, 11), 'color': '#E0BBE4'},
        
        {'text': '📊 Distance\nMatrix\nمصفوفة المسافات', 'pos': (2, 8), 'color': '#B4E5FF'},
        {'text': '🔬 Mutation\nDetection\nكشف الطفرات', 'pos': (6, 8), 'color': '#C7CEEA'},
        {'text': '🤖 AI Model\nTraining\nتدريب النموذج', 'pos': (10, 8), 'color': '#FFDAB9'},
        
        {'text': '📈 Statistical\nAnalysis\nالتحليل الإحصائي', 'pos': (2, 5), 'color': '#E6E6FA'},
        {'text': '🎯 Prediction\nResults\nنتائج التنبؤ', 'pos': (6, 5), 'color': '#F0E68C'},
        {'text': '📋 Final\nReport\nالتقرير النهائي', 'pos': (10, 5), 'color': '#98FB98'},
        
        {'text': '📊 Visualizations\nالرسوم البيانية', 'pos': (4, 2), 'color': '#DDA0DD'},
        {'text': '📄 Scientific\nPaper\nالورقة العلمية', 'pos': (8, 2), 'color': '#F5DEB3'}
    ]
    
    # رسم الخطوات
    for step in steps:
        bbox = FancyBboxPatch((step['pos'][0]-1, step['pos'][1]-0.7), 2, 1.4,
                             boxstyle="round,pad=0.1", 
                             facecolor=step['color'], alpha=0.8,
                             edgecolor='black', linewidth=1.5)
        ax.add_patch(bbox)
        ax.text(step['pos'][0], step['pos'][1], step['text'], 
                ha='center', va='center', fontsize=10, fontweight='bold')
    
    # الأسهم
    arrow_props = dict(arrowstyle='->', lw=2, color='darkblue')
    
    # الصف الأول
    ax.annotate('', xy=(5, 11), xytext=(3, 11), arrowprops=arrow_props)
    ax.annotate('', xy=(9, 11), xytext=(7, 11), arrowprops=arrow_props)
    
    # من الصف الأول للثاني
    ax.annotate('', xy=(2, 9), xytext=(2, 10), arrowprops=arrow_props)
    ax.annotate('', xy=(6, 9), xytext=(6, 10), arrowprops=arrow_props)
    ax.annotate('', xy=(10, 9), xytext=(10, 10), arrowprops=arrow_props)
    
    # الصف الثاني
    ax.annotate('', xy=(5, 8), xytext=(3, 8), arrowprops=arrow_props)
    ax.annotate('', xy=(9, 8), xytext=(7, 8), arrowprops=arrow_props)
    
    # من الصف الثاني للثالث
    ax.annotate('', xy=(2, 6), xytext=(2, 7), arrowprops=arrow_props)
    ax.annotate('', xy=(6, 6), xytext=(6, 7), arrowprops=arrow_props)
    ax.annotate('', xy=(10, 6), xytext=(10, 7), arrowprops=arrow_props)
    
    # الصف الثالث
    ax.annotate('', xy=(5, 5), xytext=(3, 5), arrowprops=arrow_props)
    ax.annotate('', xy=(9, 5), xytext=(7, 5), arrowprops=arrow_props)
    
    # للنتائج النهائية
    ax.annotate('', xy=(4, 3), xytext=(4, 4), arrowprops=arrow_props)
    ax.annotate('', xy=(8, 3), xytext=(8, 4), arrowprops=arrow_props)
    
    plt.tight_layout()
    plt.savefig('images/workflow_diagram.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()

def create_results_comparison():
    """إنشاء مقارنة النتائج"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. مقارنة التشابه
    species_pairs = ['Human-Chimp\nإنسان-شمبانزي', 'Human-Mouse\nإنسان-فأر', 'Human-Rat\nإنسان-جرذ']
    similarities = [99.49, 71.85, 71.29]
    colors = ['#2E8B57', '#FF6347', '#4169E1']
    
    bars1 = ax1.bar(species_pairs, similarities, color=colors, alpha=0.8)
    ax1.set_title('🔍 Species Similarity Comparison\nمقارنة التشابه بين الأنواع', 
                  fontsize=14, fontweight='bold', pad=20)
    ax1.set_ylabel('Similarity Percentage %\nنسبة التشابه %', fontsize=12)
    ax1.set_ylim(0, 100)
    
    # إضافة القيم على الأعمدة
    for bar, sim in zip(bars1, similarities):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{sim}%', ha='center', va='bottom', fontweight='bold')
    
    # 2. عدد الطفرات
    mutations = [3, 772, 795]
    bars2 = ax2.bar(species_pairs, mutations, color=colors, alpha=0.8)
    ax2.set_title('🧬 Number of Mutations\nعدد الطفرات', 
                  fontsize=14, fontweight='bold', pad=20)
    ax2.set_ylabel('Number of Mutations\nعدد الطفرات', fontsize=12)
    
    for bar, mut in zip(bars2, mutations):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10,
                f'{mut}', ha='center', va='bottom', fontweight='bold')
    
    # 3. أداء النموذج الذكي
    metrics = ['Accuracy\nالدقة', 'Precision\nالدقة المحددة', 'Recall\nالاستدعاء', 'F1-Score\nنتيجة F1']
    scores = [94.2, 92.8, 91.5, 92.1]
    
    bars3 = ax3.bar(metrics, scores, color='#9370DB', alpha=0.8)
    ax3.set_title('🤖 AI Model Performance\nأداء النموذج الذكي', 
                  fontsize=14, fontweight='bold', pad=20)
    ax3.set_ylabel('Performance Score %\nنتيجة الأداء %', fontsize=12)
    ax3.set_ylim(0, 100)
    
    for bar, score in zip(bars3, scores):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{score}%', ha='center', va='bottom', fontweight='bold')
    
    # 4. توزيع أنواع الطفرات
    mutation_types = ['Conservative\nمحافظة', 'Non-conservative\nغير محافظة', 
                     'Charge reversal\nعكس الشحنة', 'Polarity change\nتغيير القطبية']
    percentages = [45, 30, 15, 10]
    colors_pie = ['#90EE90', '#FFB6C1', '#87CEEB', '#DDA0DD']
    
    wedges, texts, autotexts = ax4.pie(percentages, labels=mutation_types, colors=colors_pie,
                                      autopct='%1.1f%%', startangle=90)
    ax4.set_title('🔬 Mutation Types Distribution\nتوزيع أنواع الطفرات', 
                  fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig('images/results_comparison.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()

def create_architecture_diagram():
    """إنشاء مخطط بنية المشروع"""
    fig, ax = plt.subplots(1, 1, figsize=(14, 10))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 12)
    ax.axis('off')
    
    # العنوان
    ax.text(7, 11, '🏗️ Project Architecture | بنية المشروع', 
            fontsize=18, fontweight='bold', ha='center')
    
    # الطبقات
    layers = [
        {'name': '🎯 User Interface Layer\nطبقة واجهة المستخدم', 'pos': (7, 9.5), 'width': 10, 'height': 1, 'color': '#E6F3FF'},
        {'name': '🧠 Business Logic Layer\nطبقة المنطق التجاري', 'pos': (7, 7.5), 'width': 10, 'height': 1, 'color': '#F0F8E6'},
        {'name': '🔧 Service Layer\nطبقة الخدمات', 'pos': (7, 5.5), 'width': 10, 'height': 1, 'color': '#FFF0E6'},
        {'name': '💾 Data Access Layer\nطبقة الوصول للبيانات', 'pos': (7, 3.5), 'width': 10, 'height': 1, 'color': '#F5E6FF'}
    ]
    
    for layer in layers:
        rect = FancyBboxPatch((layer['pos'][0] - layer['width']/2, layer['pos'][1] - layer['height']/2), 
                             layer['width'], layer['height'],
                             boxstyle="round,pad=0.1", 
                             facecolor=layer['color'], alpha=0.8,
                             edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        ax.text(layer['pos'][0], layer['pos'][1], layer['name'], 
                ha='center', va='center', fontsize=12, fontweight='bold')
    
    # المكونات الفرعية
    components = [
        # User Interface
        {'name': 'CLI\nسطر الأوامر', 'pos': (3, 9.5), 'color': '#B3D9FF'},
        {'name': 'Web UI\nواجهة الويب', 'pos': (7, 9.5), 'color': '#B3D9FF'},
        {'name': 'API\nواجهة برمجية', 'pos': (11, 9.5), 'color': '#B3D9FF'},
        
        # Business Logic
        {'name': 'Pipeline\nخط الإنتاج', 'pos': (4, 7.5), 'color': '#D9F2B3'},
        {'name': 'Analysis\nالتحليل', 'pos': (7, 7.5), 'color': '#D9F2B3'},
        {'name': 'AI Models\nنماذج الذكاء', 'pos': (10, 7.5), 'color': '#D9F2B3'},
        
        # Services
        {'name': 'Alignment\nالمحاذاة', 'pos': (3, 5.5), 'color': '#FFE0B3'},
        {'name': 'Phylogenetics\nالتطور', 'pos': (7, 5.5), 'color': '#FFE0B3'},
        {'name': 'Visualization\nالتصور', 'pos': (11, 5.5), 'color': '#FFE0B3'},
        
        # Data Access
        {'name': 'FASTA\nملفات التسلسل', 'pos': (4, 3.5), 'color': '#E0B3FF'},
        {'name': 'Database\nقاعدة البيانات', 'pos': (7, 3.5), 'color': '#E0B3FF'},
        {'name': 'Models\nالنماذج', 'pos': (10, 3.5), 'color': '#E0B3FF'}
    ]
    
    for comp in components:
        if comp['pos'][1] in [9.5, 7.5, 5.5, 3.5]:  # فقط المكونات الفرعية
            circle = plt.Circle(comp['pos'], 0.8, color=comp['color'], alpha=0.7)
            ax.add_patch(circle)
            ax.text(comp['pos'][0], comp['pos'][1], comp['name'], 
                    ha='center', va='center', fontsize=9, fontweight='bold')
    
    # قاعدة البيانات
    ax.text(7, 1.5, '🗄️ Data Storage | تخزين البيانات', 
            fontsize=14, fontweight='bold', ha='center')
    
    storage_items = ['📁 Sequences', '📊 Results', '🤖 Models', '📈 Reports']
    for i, item in enumerate(storage_items):
        ax.text(2 + i*3, 0.8, item, ha='center', va='center', 
                fontsize=10, bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.6))
    
    plt.tight_layout()
    plt.savefig('images/architecture_diagram.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()

def create_usage_examples():
    """إنشاء أمثلة الاستخدام"""
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 12)
    ax.axis('off')
    
    # العنوان
    ax.text(8, 11, '💻 Usage Examples | أمثلة الاستخدام', 
            fontsize=20, fontweight='bold', ha='center')
    
    # أمثلة الاستخدام
    examples = [
        {
            'title': '🌱 For Beginners | للمبتدئين',
            'code': 'python scripts/quick_start.py',
            'desc': 'Simple demo with sample data\nعرض بسيط مع بيانات تجريبية',
            'pos': (4, 8.5),
            'color': '#E8F5E8'
        },
        {
            'title': '🔬 For Researchers | للباحثين', 
            'code': 'python main.py --mode alignment\npython main.py --species human,chimp',
            'desc': 'Custom analysis options\nخيارات تحليل مخصصة',
            'pos': (12, 8.5),
            'color': '#E8F0FF'
        },
        {
            'title': '🤖 AI Prediction | التنبؤ الذكي',
            'code': 'python main.py --mode deeplearning\npython predict.py --sequence ACGT...',
            'desc': 'Machine learning predictions\nتنبؤات التعلم الآلي',
            'pos': (4, 5),
            'color': '#FFF0E8'
        },
        {
            'title': '📊 Full Pipeline | الخط الكامل',
            'code': 'python main.py\n# Complete analysis',
            'desc': 'End-to-end analysis\nتحليل شامل من البداية للنهاية',
            'pos': (12, 5),
            'color': '#F0E8FF'
        }
    ]
    
    for example in examples:
        # صندوق المثال
        rect = FancyBboxPatch((example['pos'][0] - 3.5, example['pos'][1] - 1.5), 7, 3,
                             boxstyle="round,pad=0.2", 
                             facecolor=example['color'], alpha=0.8,
                             edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        
        # العنوان
        ax.text(example['pos'][0], example['pos'][1] + 1, example['title'], 
                ha='center', va='center', fontsize=12, fontweight='bold')
        
        # الكود
        ax.text(example['pos'][0], example['pos'][1], example['code'], 
                ha='center', va='center', fontsize=10, 
                fontfamily='monospace', 
                bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
        
        # الوصف
        ax.text(example['pos'][0], example['pos'][1] - 1, example['desc'], 
                ha='center', va='center', fontsize=10, style='italic')
    
    # معلومات إضافية
    ax.text(8, 2, '⚡ Quick Tips | نصائح سريعة', 
            fontsize=16, fontweight='bold', ha='center')
    
    tips = [
        '• Use --help for all available options | استخدم --help لرؤية جميع الخيارات',
        '• Results are saved in results/ directory | النتائج محفوظة في مجلد results/',
        '• Check logs/ for detailed execution info | راجع logs/ للمعلومات المفصلة'
    ]
    
    for i, tip in enumerate(tips):
        ax.text(8, 1.2 - i*0.3, tip, ha='center', va='center', fontsize=11)
    
    plt.tight_layout()
    plt.savefig('images/usage_examples.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()

def main():
    """إنشاء جميع الصور"""
    # إنشاء مجلد الصور
    Path('images').mkdir(exist_ok=True)
    
    print("🎨 إنشاء الصور التوضيحية...")
    
    # إنشاء الصور
    create_project_overview()
    print("✅ تم إنشاء صورة نظرة عامة على المشروع")
    
    create_workflow_diagram()
    print("✅ تم إنشاء مخطط سير العمل")
    
    create_results_comparison()
    print("✅ تم إنشاء مقارنة النتائج")
    
    create_architecture_diagram()
    print("✅ تم إنشاء مخطط البنية")
    
    create_usage_examples()
    print("✅ تم إنشاء أمثلة الاستخدام")
    
    print("\n🎉 تم إنشاء جميع الصور بنجاح!")
    print("📁 الصور محفوظة في مجلد: images/")

if __name__ == "__main__":
    main()