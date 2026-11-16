"""Run the Sage source with consistent argv/__name__ across Sage 9.5 and 10.8."""
import sys
from pathlib import Path
from sage.repl.preparse import preparse

source = Path(__file__).with_name('verify_and_graphs.sage')
sys.argv[0] = str(source)
exec(compile(preparse(source.read_text()), str(source), 'exec'),
     {'__name__': '__main__', '__file__': str(source)})
