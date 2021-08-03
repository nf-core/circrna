'use strict';
const spawn = require('buffered-spawn');
const fail = require('./git-author').fail;

/**
 * 运行 git blame 命令
 *
 * @param   {Vinyl}   file      要检查blame信息的文件
 * @param   {Boolean} showEmail 是否使用后email代替用户名
 * @returns {Promise} 命令运行结果，Buffer
 */
function runBlame (file) {
	const args = [
		'--no-pager',
		'blame',
		'-w',
		'-C',
		'-M',
		'-p',
		'--contents',
		'-',
	];
	let write;

	if (file.isStream()) {
		write = (stdin) => {
			file.contents.pipe(stdin);
		};
	} else if (file.isBuffer()) {
		write = (stdin) => {
			stdin.write(file.contents);
			stdin.end();
		};
	} else {
		args.length = 6;
	}

	args.push('--', file.path);

	const blame = spawn('git', args, {
		crossSpawn: false,
		cwd: file.cwd,
	});

	if (write) {
		write(blame.cp.stdin);
		blame.cp.stdin.on('error', () => {
			// console.error
		});
	}

	return blame;
}

function parsePerson (data, result) {
	data.replace(/^\w+-(\w+) (.*)$/gm, (s, key, value) => {
		if (key === 'mail') {
			value = value.replace(/^<(.*)>$/, '$1');
		} else if (key === 'time') {
			value = value - 0;
		}
		result[key] = value;
		return '';
	});
	return result;
}

function parseLine (props, result) {
	props.replace(/^(\w+) (.*)(\n(?:\1-\w+ .*\n)+)/gm, (s, role, name, props) => {
		result.rev[role] = parsePerson(props, {
			name,
		});
		return '';
	}).replace(/^summary (.*)$/igm, (s, summary) => {
		result.rev.summary = summary;
		return '';
	}).replace(/^previous (\w+) (.*)$/igm, (s, hash, filename) => {
		result.previous = {
			hash,
			filename,
		};
		return '';
	}).replace(/^(\S+) (.*)$/gm, (s, key, value) => {
		result[key] = value;
		return '';
	});
	return result;
}
/**
 * 将git blame 命令运行结果转换为行号作为下标的数组
 *
 * @param {string}   data git blame命令返回的buffer信息
 * @returns {Blame[]}     数组，里面是每行代码的作者信息对象
 */
function parseBlame (data) {
	const revCache = {};
	const result = [];
	data.replace(/^(\w{40,}) (\d+) (\d+)(?: \d+)*\n((?:\S*.*\n)*?)\t(.*)$/gm, (s, hash, originalLine, finalLine, props, content) => {
		originalLine = originalLine - 0;
		finalLine = finalLine - 0;
		result[finalLine] = parseLine(props, {
			originalLine,
			finalLine,
			content,
			rev: revCache[hash] || (revCache[hash] = {
				hash,
			}),
		});
		return '';
	});
	return result;
}

/**
 * 将git blame 命令运行结果转换为行号作为下标的数组
 *
 * @param {Buffer?}   data git blame命令返回的buffer信息
 * @returns {Blame[]}      数组，里面是每行代码的作者信息对象
 */
function parseReaslt (data) {
	data = data && data.toString();
	if (data) {
		return parseBlame(data);
	}
}

/**
 * 运行 git blame 命令
 *
 * @async
 * @param   {Vinyl}   file      要检查blame信息的文件
 * @returns {Promise<Blame[]>} 命令运行结果，数组，里面是每行代码的作者信息对象
 */
function getFileBlame (file) {
	if (fail[file.cwd]) {
		return Promise.reject(fail[file.cwd]);
	}
	return runBlame(file).then(cp => parseReaslt(cp.stdout));
}

module.exports = getFileBlame;
