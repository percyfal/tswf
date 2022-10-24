from tswf.cli import cli


def test_cli(runner):
    result = runner.invoke(cli, [])
    print(result.output)
