# Хеш-функция Стрибог (ГОСТ 34.11-2012)

[![ГОСТ 34.11-2012](https://img.shields.io/badge/ГОСТ-34.11--2012-blue?style=plastic)](https://meganorm.ru/Data2/1/4293788/4293788459.pdf)
![C++20](https://img.shields.io/badge/std-C++20-purple?style=plastic)
![Google Code Style](https://img.shields.io/badge/style-Google-purple?style=plastic)
[![Benchmark](https://img.shields.io/badge/benchmark-hashing%20race-brightgreen?style=plastic)](doc/benchmarks.md)
[![GitHub license](https://img.shields.io/github/license/gdaneek/streebog-hash?style=plastic)](https://github.com/gdaneek/streebog-hash/blob/master/LICENSE)

Высокопроизводительная реализация хеш-функции **Стрибог** согласно российским криптографическим стандартам ГОСТ 34.11-2012 (ГОСТ 34.11-2018).

---

## ⚙️ Установка

1. Клонируйте репозиторий:
```bash
git clone https://github.com/gdaneek/streebog-hash.git
cd streebog-hash
```

2. Сборка проекта с помощью CMake:
```bash
mkdir build
cd build
cmake ..
cmake --build .
```

## 📦 Результаты сборки

После сборки вы получите следующие компоненты:

- 🔹 `stbg512` — утилита для вычисления 512-битного хеша
- 🔹 `stbg256` — утилита для вычисления 256-битного хеша
- 🔹 `stbg` — универсальная утилита для обоих режимов (256/512)
- 📚 `libstreebog.a` — статическая библиотека
- 🧪 Тесты — основаны на официальных примерах из стандарта

## 🧑‍💻 Документация разработчика

Документация находится в папке [`doc/code`](doc/code) или может быть сгенерирована с помощью **Doxygen**:

```bash
cd doc
doxygen

```

## 📊 Бенчмарки

См. файл [`doc/benchmarks.md`](doc/benchmarks.md) для подробных результатов производительности.

## 🔧 Продвинутые настройки

- ✅ Для упрощения работы с функцией как с STL-контейнером доступны специальные обёртки. Включаются с помощью макроса:

```cpp
#define STREEBOG_ENABLE_WRAPPERS
```

---

## 📄 Лицензия

Проект распространяется под лицензией **GNU GPL**. Подробнее см. в [LICENSE](LICENSE).
