from asyncio import new_event_loop
from os.path import join
from tempfile import TemporaryDirectory

loop = browser = None

async def render(s, t, width, height, scale):
    page = await browser.newPage()
    await page.setViewport({'deviceScaleFactor': scale, 'width': width, 'height': height})
    await page.goto(f'file://{s}')
    element = await page.querySelector('svg')
    await element.screenshot({'path': t})
    await page.close()

def svg2png(svg: str, width: int = 1000, height: int = 1000, scale: float = 10.):
    global loop, browser

    if loop is None:  # lazy browser launcher
        from pyppeteer import launch

        loop = new_event_loop()
        browser = loop.run_until_complete(launch())
    elif browser is None:
        raise ImportError('pyppeteer initialization failed')

    with TemporaryDirectory() as tmpdir:
        with open(s := join(tmpdir, 'input.svg'), 'w') as f:
            f.write(svg)

        loop.run_until_complete(render(s, (t := join(tmpdir, 'output.png')), width, height, scale))

        with open(t, 'rb') as f:
            return f.read()

__all__ = ['svg2png']
