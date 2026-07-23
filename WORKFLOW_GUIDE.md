# GitHub 업로드 가이드 (폴더를 git 저장소로 만들기 → GitHub 연결 → 커밋/push)

이 문서는 `pipeline_synthesis` 폴더를 git 저장소로 만들고 GitHub과 연결해서 커밋/push하기까지의
과정을 정리한 것이다. (문서 내용을 어떻게 정리할지에 대한 방법론은 다루지 않는다 — 여기서는
순수하게 "폴더 → git → GitHub" 과정만 다룬다.)

---

## 0-1. 브랜치(branch) 구조

지금 이 저장소는 `main` 브랜치 하나만 쓰고 있고(`git branch -a`로 확인), 커밋이 일자로 쭉
이어지는 단순한 구조다. 혼자 쓰는 문서 저장소라면 **이 상태(브랜치 하나, 계속 main에 커밋)가
정상이며 가장 단순하다** — 여러 브랜치를 나누는 건 "동시에 여러 작업을 병행해야 할 때"나
"아직 확정 안 된 내용을 완성본과 분리해두고 싶을 때" 필요한 것이지, 항상 써야 하는 게 아니다.

브랜치가 실제로 유용한 상황 하나: **아직 검증이 덜 끝난 내용을 작업할 때**. 이전에 Rubus
프로젝트를 정리할 때, main에 바로 커밋했다가 "확인 필요 항목이 너무 많아서 이대로는 GitHub에
올리기 싫다"는 상황이 생겨서, 태그로 백업하고 `git reset --hard`로 되돌리는 다소 번거로운
작업을 했다(아래 참고 섹션). **애초에 draft용 브랜치에서 작업했으면 이 번거로움이 없었다.**

### 브랜치 기본 명령어

```bash
git branch                      # 현재 브랜치 목록 보기 (로컬)
git branch -a                   # 원격 브랜치까지 포함해서 보기
git branch <새브랜치이름>        # 새 브랜치 생성(전환은 안 됨)
git checkout -b <새브랜치이름>   # 새 브랜치 생성 + 바로 전환
git checkout main                # main으로 다시 전환
git branch -d <브랜치이름>       # 다 쓴 브랜치 삭제(로컬)
```

### 앞으로 쓸 수 있는 패턴: "draft 브랜치에서 검증, 끝나면 main으로"

```bash
# 1. 새 프로젝트 조사/정리를 시작할 때 draft 브랜치를 만든다
git checkout -b draft/새프로젝트이름

# 2. 이 브랜치 위에서 마음껏 커밋한다(확인 필요 항목이 있어도 상관없음 — 아직 main이 아니므로)
git add ...
git commit -m "..."

# 3. 검증이 끝나고 "이제 GitHub에 올려도 되겠다" 싶으면 main으로 합친다
git checkout main
git merge draft/새프로젝트이름     # draft의 커밋들을 main에 합침

# 4. main만 push (draft 브랜치는 로컬에만 있어도 되고, push 안 하면 GitHub엔 안 보임)
git push https://계정:토큰@github.com/.../저장소.git main:main

# 5. 다 쓴 draft 브랜치는 삭제(선택)
git branch -d draft/새프로젝트이름
```

이 방식의 장점: `git reset --hard`처럼 브랜치 히스토리를 되돌리는 위험한 작업을 할 필요가
없다. draft 브랜치는 그냥 "아직 main에 합치지 않은 별도의 작업 공간"일 뿐이라, 검증 중에
커밋을 몇 개를 하든 main에는 전혀 영향이 없고, 확신이 설 때만 `merge`로 반영하면 된다.

> 주의: `merge`는 draft 브랜치의 **모든 커밋**을 main으로 가져온다. draft에서 커밋을 여러 개
> 나눠서 했더라도(예: "1차 조사" → "확인 필요 제거" → "재작성"), merge하면 그 전체 히스토리가
> main에 들어간다. 최종 결과물 하나만 깔끔하게 main에 남기고 싶다면 `git merge --squash
> draft/이름` 후 `git commit`을 쓰면, draft의 여러 커밋이 하나의 커밋으로 합쳐져서 main에
> 들어간다(이 방법이 지금까지 해온 "확인된 내용만 새 커밋 하나로" 패턴과 가장 잘 맞는다).

---

## 0-2. 이미 되어 있는 것 (참고용)

`pipeline_synthesis` 폴더는 이미 아래 상태로 설정되어 있다:

```bash
cd /pbi-acc1/hsPark/pipeline_synthesis
git status          # 이미 git 저장소임
git remote -v       # origin이 이미 GitHub과 연결되어 있음
```

새 폴더를 처음부터 설정할 때는 아래 1~3번 과정을 그대로 따라가면 된다.

---

## 1. 폴더를 git 저장소로 만들기

```bash
cd /경로/정리할_폴더
git init
```

이 저장소 전용 커밋 계정을 설정한다(`--global` 없이 — 이 폴더에만 적용됨):

```bash
git config user.name "이름"
git config user.email "이메일주소"
```

(글로벌 설정이 이미 되어 있다면 이 단계는 생략 가능. `pipeline_synthesis`는 이 저장소에만
`user.name "hsPark"` / `user.email "pbicnu@gmail.com"`로 로컬 설정되어 있다.)

---

## 2. GitHub에 원격 저장소 만들고 연결하기

### 2-1. GitHub 웹사이트에서 새 저장소 생성

- github.com → New repository
- Private 선택(민감한 내용이면), README 없이 **빈 저장소로 생성** (로컬에 이미 커밋이 있을
  경우 README를 자동 생성하면 충돌 남)
- 생성 후 나오는 저장소 URL을 기록: `https://github.com/<계정>/<저장소이름>.git`

### 2-2. 로컬 저장소에 원격 연결

```bash
git remote add origin https://github.com/<계정>/<저장소이름>.git
git branch -M main
git remote -v   # 확인
```

---

## 3. 커밋

```bash
git add <파일들>
git commit -m "커밋 메시지"
```

---

## 4. GitHub으로 push (인증)

이 서버는 브라우저 로그인이 불가능하고 `gh` CLI도 없어서, **Personal Access Token(PAT)을
push 커맨드에 직접 넣는 방식**을 쓴다.

### 4-1. Fine-grained PAT 발급

GitHub → Settings → Developer settings → Personal access tokens → Fine-grained tokens:

- **Repository access**: 해당 저장소만 선택
- **Permissions → Repository permissions → Contents**: **Read and write**로 변경
  (이것만 바꾸면 push 가능. 나머지 권한은 기본값 그대로 둬도 됨)
- 만료기간은 짧게(예: 7일) 설정 권장

### 4-2. push

처음 push할 때(원격 브랜치가 아직 없을 때):

```bash
git push -u https://<GitHub계정>:<발급받은_토큰>@github.com/<계정>/<저장소이름>.git main
```

이후 push할 때:

```bash
git push https://<GitHub계정>:<발급받은_토큰>@github.com/<계정>/<저장소이름>.git main:main
```

### 4-3. push 후 정리

- push가 끝나면 **바로 해당 토큰을 GitHub에서 revoke(삭제)**한다 — 명령어 자체에 토큰이
  그대로 노출되므로.
- 원한다면 remote URL에서 토큰 흔적을 지우기 위해 다시 깨끗한 URL로 재설정할 수 있다(토큰이
  git config에 저장되는 건 아니지만, 습관적으로 정리하고 싶다면):
  ```bash
  git remote set-url origin https://github.com/<계정>/<저장소이름>.git
  ```

---

## 5. 이후 작업 (파일 수정 → 커밋 → push 반복)

한 번 1~2번 설정이 끝나면, 그 다음부터는 3~4번만 반복하면 된다:

```bash
git add <바뀐 파일들>
git commit -m "메시지"
git push https://<계정>:<새로_발급한_토큰>@github.com/<계정>/<저장소이름>.git main:main
```

---

## 참고: 특정 커밋을 GitHub에는 올리지 않고 로컬에만 남기고 싶을 때

이미 로컬에 커밋해둔 내용 중 일부를 GitHub에는 올리지 않고 싶다면(예: 아직 검증이 덜 된
내용), git push는 브랜치의 커밋 히스토리를 통째로 보낸다는 점을 유의해야 한다 — 그 커밋
위에 새 커밋을 얹어서 push하면 이전 커밋도 같이 올라간다.

로컬에는 남기고 GitHub에는 안 올리려면:

```bash
# 1. 문제의 커밋을 태그로 백업(로컬에만 존재, push 안 하면 GitHub에 안 올라감)
git tag backup-이름 <해당커밋해시>

# 2. 브랜치를 그 이전(마지막으로 이미 push된) 커밋으로 되돌림
git reset --hard <마지막으로_push한_커밋해시>

# 3. 필요하면 내용을 다시 정리해서 새 커밋
# 4. push (이번엔 백업 태그로 남긴 커밋은 히스토리에 없으므로 GitHub에 안 올라감)
```

백업 태그는 나중에 이렇게 다시 확인할 수 있다:

```bash
git tag -l                      # 백업 태그 목록
git show backup-이름:파일명      # 백업 시점 특정 파일 내용 보기
```
