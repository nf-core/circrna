'use strict';

const bufferedSpawn = require('../index');
const expect = require('expect.js');

const isWin = process.platform === 'win32';

describe('buffered-spawn', () => {
    describe('api', () => {
        it('should handle optional args & options', (next) => {
            bufferedSpawn('echo', (err) => {
                expect(err).to.not.be.ok();
                next();
            });
        });

        it('should handle optional args', (next) => {
            bufferedSpawn(`${__dirname}/fixtures/hello`, {
                stdio: ['pipe', 'ignore', 'ignore'],
            }, (err, stdout) => {
                expect(err).to.not.be.ok();
                expect(stdout).to.be('');

                bufferedSpawn(`${__dirname}/fixtures/hello`, null, {
                    stdio: ['pipe', 'ignore', 'ignore'],
                }, (err, stdout) => {
                    expect(err).to.not.be.ok();
                    expect(stdout).to.be('');

                    next();
                });
            });
        });

        it('should handle optional options', (next) => {
            bufferedSpawn('node', [
                `${__dirname}/fixtures/echo`,
                'foo',
            ], (err, stdout) => {
                expect(err).to.not.be.ok();
                expect(stdout.trim()).to.equal('foo');

                bufferedSpawn('node', [
                    `${__dirname}/fixtures/echo`,
                    'foo',
                ], null, (err, stdout) => {
                    expect(err).to.not.be.ok();
                    expect(stdout.trim()).to.equal('foo');

                    next();
                });
            });
        });

        it('should pass arguments to node\'s spawn', (next) => {
            bufferedSpawn('node', ['simple'], { cwd: `${__dirname}/fixtures` }, (err, stdout, stderr) => {
                expect(err).to.not.be.ok();
                expect(stdout).to.equal('i am being printed on stdout');
                expect(stderr).to.equal('i am being printed on stderr');

                next();
            });
        });

        it('should allow node\'s spawn\'s stdout to be ignored & inherited', (next) => {
            bufferedSpawn('node', ['simple'], {
                cwd: `${__dirname}/fixtures`,
                stdio: ['pipe', 'ignore', 'pipe'],
            }, (err, stdout, stderr) => {
                expect(err).to.not.be.ok();
                expect(stdout).to.equal('');
                expect(stderr).to.equal('i am being printed on stderr');

                bufferedSpawn('node', ['simple'], {
                    cwd: `${__dirname}/fixtures`,
                    stdio: ['pipe', 1, 'pipe'],
                }, (err, stdout, stderr) => {
                    expect(err).to.not.be.ok();
                    expect(stdout).to.equal('');
                    expect(stderr).to.equal('i am being printed on stderr');

                    next();
                });
            });
        });

        it('should allow node\'s spawn\'s stderr to be ignored & inherited', (next) => {
            bufferedSpawn('node', ['simple'], {
                cwd: `${__dirname}/fixtures`,
                stdio: ['pipe', 'pipe', 'ignore'],
            }, (err, stdout, stderr) => {
                expect(err).to.not.be.ok();
                expect(stdout).to.equal('i am being printed on stdout');
                expect(stderr).to.equal('');

                bufferedSpawn('node', ['simple'], {
                    cwd: `${__dirname}/fixtures`,
                    stdio: ['pipe', 'pipe', 2],
                }, (err, stdout, stderr) => {
                    expect(err).to.not.be.ok();
                    expect(stdout).to.equal('i am being printed on stdout');
                    expect(stderr).to.equal('');

                    next();
                });
            });
        });

        it('should work with promises', () => {
            return bufferedSpawn('node', [
                `${__dirname}/fixtures/echo`,
                'foo',
            ])
            .then((output) => {
                expect(output.stdout.trim()).to.equal('foo');
                expect(output.stderr.trim()).to.equal('');

                return bufferedSpawn('node', [`${__dirname}/fixtures/fail`]);
            })
            .then(() => {
                throw new Error('Should have failed');
            }, (err) => {
                expect(err).to.be.an(Error);
                expect(err.status).to.equal(25);
                expect(err.stdout).to.equal('stdout fail');
                expect(err.stderr).to.equal('stderr fail');
                expect(err.details).to.equal(err.stderr);
            });
        });

        it('should give access to the underlying child process', () => {
            const cp = bufferedSpawn('echo', () => {});

            expect(cp.kill).to.be.a('function');

            const promise = bufferedSpawn('echo');

            expect(promise.cp.kill).to.be.a('function');
        });
    });

    it('should buffer stdout & stderr', () => {
        return bufferedSpawn('node', [`${__dirname}/fixtures/simple`])
        .then((output) => {
            expect(output.stdout).to.equal('i am being printed on stdout');
            expect(output.stderr).to.equal('i am being printed on stderr');
        });
    });

    it('should expand using PATH_EXT properly', () => {
        if (!isWin) {
            return Promise.resolve();
        }

        return bufferedSpawn(`${__dirname}/fixtures/foo`)  // Should expand to foo.bat
        .then((output) => {
            expect(output.stdout.trim()).to.equal('foo');

            return bufferedSpawn(`${__dirname}/fixtures/foo`, [], { crossSpawn: false });
        })
        .then(() => {
            throw new Error('Should have failed');
        }, (err) => {
            expect(err).to.be.an(Error);
            expect(err.code).to.be('ENOENT');
        });
    });

    it('should handle multibyte properly', () => {
        return bufferedSpawn('node', [`${__dirname}/fixtures/multibyte`])
        .then((output) => {
            expect(output.stdout).to.equal('こんにちは');
            expect(output.stderr).to.equal('こんにちは');
        });
    });

    it('should fail on error code != 0 and still give stdout/stderr', () => {
        return bufferedSpawn('node', [`${__dirname}/fixtures/fail`])
        .then(() => {
            throw new Error('Should have failed');
        }, (err) => {
            expect(err).to.be.an(Error);
            expect(err.status).to.equal(25);
            expect(err.stdout).to.equal('stdout fail');
            expect(err.stderr).to.equal('stderr fail');
            expect(err.details).to.equal(err.stderr);
        });
    });
});
