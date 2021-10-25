from urllib.parse import urlparse
from bs4 import BeautifulSoup 
import requests, os

# URL
URL = "https://en.wiktionary.org"

page = requests.get(URL + "/wiki/Category:Japanese_appendices")
soup = BeautifulSoup(page.content, "html.parser")

# <div> containing JLPT links 
results = soup.find_all("div", class_="mw-category-group")[0]

# get links to each JLPT page 
results = [link.get("href") for link in results.find_all("a")[:5]]

# compile links into full (not relative) URLs into dict 
JLPT_URLs = {link[-2:] : urlparse(URL + link).geturl() for link in results}

print("Parsed the following urls: ")
for k, v in JLPT_URLs.items():
    print(k) 
    print(v + " \n")

def parse_jlpt(key: str, dest: str, card_suffix="#card-reverse", show_test=True):
    """Parse vocabulary lists in link corresponding to `key` of `JLPT_URLs`, and save to a file at `dest`

    Args:
        key (str): one of 'N1', 'N2', ... 'N5'
        dest (str): destination of output file
        card_suffix (optional, str): suffix to add to end of each Japanese vocabulary. Default is "card-reverse", which creates Front-Back cards when using Obsidian and the Flashcards plugin.
        show_test (optional, bool): whether to print first 6 lines of output file. Defaults to True.
    """
    
    if not os.path.isdir(dest):
        raise ValueError("Output destination %s is not a valid directory" % dest)
    if not isinstance(card_suffix, str):
        raise ValueError(
            "`card_suffix` must be a string, not <{0}>".format(type(card_suffix))
        )
        
    out_path = dest + key + "_vocab.md"
    if os.path.isfile(out_path):
        print("File already exists at < %s >. Skipping." % out_path)
        return None 
    
    url = JLPT_URLs[key]
    print("Reading... ", url)
        
    page = requests.get(url)
    soup = BeautifulSoup(page.content, "html.parser")
    
    results = soup.find_all("div", class_="mw-parser-output")[0]
    
    # vocabualry can be in unordered <ul> or ordered <ol> lists
    ord_ = results.find_all("ol")
    un_ = results.find_all("ul")
    
    try: 
        assert len(ord_) > 1 or len(un_) > 1
    except:
        print("No vocabulary found using <ol> and <ul> tags: \n %s" % url)
        
    if len(ord_) > 0:
        vocab = ord_ 
    else:
        vocab = un_ 
        
    if len(vocab) > 1:
        for v in vocab:            
            if "</span> -" in str(v.find("li")): 
                vocab = v 
                break 
    else: 
        vocab = vocab[0]
    
    vocab = [s.text for s in vocab.find_all("li")]
    
    out_file = open(out_path, "w", encoding="utf8")
    
    line_fmt = "{0} {1}\n{2}\n\n"
    for line in vocab:
        if " -" in line:
            v = line.split(" -")
            out_file.write(line_fmt.format(v[0], card_suffix, v[1]))
        else: continue 
        
    out_file.close()
    
    print("Successfully parsed < %s > vocabulary and saved to < %s >" % (key, out_path))
    
    if show_test:
        print("First 3 lines: ")
        with open(out_path, "r", encoding="utf8") as f:
            lines = f.readlines()[:9]
            for i in range(0, 9, 3):
                print(lines[i:i+2])


out_path = r"C:/Users/delbe/Downloads/wut/wut/Post_grad/UBC/Research/records/obsidian notes/Japanese/"    

for key in JLPT_URLs.keys():
    parse_jlpt(key, out_path)
